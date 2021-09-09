module task_manager_module
  !> A task manager determines the distribution of tasks to workers with
  !> respect to the tasks' demands.

  use libatoms_module, only : initialise, finalise, print_title, print, &
      operator(//), system_abort, inoutput, optional_default, OUTPUT
  use linearalgebra_module, only : heap_sort
  use MPI_context_module, only : MPI_Context, print
  use ScaLAPACK_module, only : Scalapack, print

  implicit none
  private

  public :: task_manager_type, task_type, worker_type
  public :: task_manager_init_tasks, task_manager_init_workers, task_manager_deallocate
  public :: task_manager_add_task, task_manager_add_worker
  public :: task_manager_distribute_tasks
  public :: task_manager_print, task_manager_export_distribution

  integer, parameter :: UNDEFINED = -1

  type task_type
    integer :: index = UNDEFINED !> initial task index
    integer :: memory = 0 !> memory needed for this task
    integer :: worker_id = UNDEFINED !> ID of worker assigned to this task
  end type task_type

  type worker_type
    integer :: index = UNDEFINED !> initial worker index
    integer :: n_tasks = 0 !> number of tasks assigned to this worker
    integer :: n_padding = 0 !> number of padding rows to reach unified workload
    integer, dimension(:), allocatable :: task_ids !> IDs of tasks assigned to this worker
  end type worker_type

  type task_manager_type
    logical :: active = .false. !> whether task manager is used
    logical :: distributed = .false. !> whether tasks have been distributed
    integer :: my_worker_id = 0 !> compute only tasks for this worker
    integer :: unified_workload = 0 !> highest workload among workers, target for padding

    integer :: n_tasks = 0 !> number of tasks
    integer :: n_workers = 0 !> number of workers

    type(MPI_Context) :: MPI_obj
    type(ScaLAPACK) :: ScaLAPACK_obj

    type(task_type), dimension(:), allocatable :: tasks !> tasks to compute
    type(worker_type), dimension(:), allocatable :: workers !> workers to compute tasks
  end type task_manager_type

  contains

  !> Add task to list
  subroutine task_manager_add_task(this, memory, worker_id)
    type(task_manager_type), intent(inout) :: this
    integer, intent(in) :: memory !> memory needed for new task
    integer, intent(in), optional :: worker_id !> preset worker for this task

    this%n_tasks = this%n_tasks + 1
    if (this%n_tasks > size(this%tasks)) call system_abort("More tasks added than allocated.")
    associate(task => this%tasks(this%n_tasks))
        task%index = this%n_tasks
        task%memory = memory
        if (present(worker_id)) task%worker_id = worker_id
    end associate
  end subroutine task_manager_add_task

  !> Add worker to list
  subroutine task_manager_add_worker(this)
    type(task_manager_type), intent(inout) :: this

    this%n_workers = this%n_workers + 1
    if (this%n_workers > size(this%workers)) call system_abort("More workers added than allocated.")
    associate(worker => this%workers(this%n_workers))
        worker%index = this%n_workers
        worker%n_tasks = 0
    end associate
  end subroutine task_manager_add_worker

  !> @brief Distribute tasks among workers for minimal memory usage
  !> Assertions: %tasks(:)%memory is set, %n_workers > 0
  subroutine task_manager_distribute_tasks(this)
    type(task_manager_type), intent(inout) :: this

    integer :: i, t, w ! indices
    integer :: n
    integer, dimension(:), allocatable :: workloads
    integer, dimension(:), allocatable :: memory_list
    integer, dimension(:), allocatable :: index_list

    allocate(workloads(this%n_workers))
    allocate(memory_list(this%n_tasks))
    allocate(index_list(this%n_tasks))

    ! copy to temp arrays, sort both by memory
    memory_list = this%tasks(1:this%n_tasks)%memory
    index_list = this%tasks(1:this%n_tasks)%index
    call heap_sort(memory_list, i_data=index_list) ! index_list sorted via memory_list

    ! count workload of preset tasks
    workloads = 0
    do t = 1, this%n_tasks
        i = index_list(t)
        w = this%tasks(i)%worker_id
        if (w /= UNDEFINED) then
          this%workers(w)%n_tasks = this%workers(w)%n_tasks + 1
          workloads(w) = workloads(w) + memory_list(t)
        end if
    end do

    ! distribute tasks to workers to minimize maximum memory usage
    do t = this%n_tasks, 1, -1
        i = index_list(t)
        if (this%tasks(i)%worker_id /= UNDEFINED) cycle ! this preset tasks
        w = minloc(workloads, 1)
        this%tasks(i)%worker_id = w
        this%workers(w)%n_tasks = this%workers(w)%n_tasks + 1
        workloads(w) = workloads(w) + memory_list(t)
    end do

    this%unified_workload = maxval(workloads)

    ! fill workers with task id lists (replication for simpler retrieval)
    do w = 1, this%n_workers
        n = this%workers(w)%n_tasks
        allocate(this%workers(w)%task_ids(n))
        this%workers(w)%n_tasks = 0
    end do
    do t = 1, this%n_tasks
        w = this%tasks(t)%worker_id
        n = this%workers(w)%n_tasks + 1
        this%workers(w)%task_ids(n) = t
        this%workers(w)%n_tasks = n
    end do

    do w = 1, this%n_workers
        this%workers(w)%n_padding = this%unified_workload - workloads(w)
    end do

    this%distributed = .true.
  end subroutine task_manager_distribute_tasks

  subroutine task_manager_show_distribution(this)
    type(task_manager_type), intent(in) :: this

    integer :: i, t, w, load, n_tasks

    do w = 1, this%n_workers
        call print("Tasks of worker "//w)
        load = 0
        n_tasks = 0
        do i = 1, this%workers(w)%n_tasks
          t = this%workers(w)%task_ids(i)
          call print(i // ": " // t // ": " // this%tasks(t)%memory)
          load = load + this%tasks(t)%memory
        end do
        call print("Summary: "//this%workers(w)%n_tasks//" tasks, "//load//" memory, "//this%workers(w)%n_padding//" padding")
    end do
  end subroutine task_manager_show_distribution

  subroutine task_manager_export_distribution(this, only_id)
    type(task_manager_type), intent(in) :: this
    integer, optional :: only_id !> only this worker writes distribution

    integer :: i, t, w
    integer :: only_id_opt
    type(Inoutput) :: file

    only_id_opt = optional_default(UNDEFINED, only_id)
    if (only_id_opt >= 0 .and. only_id_opt /= this%my_worker_id) return

    call initialise(file, "tm.dist", OUTPUT)
    call print("worker_id task_id_local task_id_global memory", file=file)
    do w = 1, this%n_workers
        do i = 1, this%workers(w)%n_tasks
          t = this%workers(w)%task_ids(i)
          call print(w // " " // i // " " // t // " " // this%tasks(t)%memory, file=file)
        end do
    end do
    call finalise(file)
  end subroutine task_manager_export_distribution

  subroutine task_manager_init_tasks(this, n_tasks)
    type(task_manager_type), intent(inout) :: this
    integer, intent(in) :: n_tasks !> number of tasks to allocate

    call print("Allocate %tasks for "//n_tasks//" tasks.")
    allocate(this%tasks(n_tasks))
  end subroutine task_manager_init_tasks

  subroutine task_manager_init_workers(this, n_workers)
    type(task_manager_type), intent(inout) :: this
    integer, intent(in) :: n_workers !> number of workers to allocate and set

    this%n_workers = n_workers
    allocate(this%workers(this%n_workers))
  end subroutine task_manager_init_workers

  subroutine task_manager_deallocate(this)
    type(task_manager_type), intent(inout) :: this

    this%distributed = .false.
    this%n_tasks = 0
    this%n_workers = 0
    if (allocated(this%tasks)) deallocate(this%tasks)
    if (allocated(this%workers)) deallocate(this%workers)
  end subroutine task_manager_deallocate

  subroutine task_manager_print(this)
    type(task_manager_type), intent(inout) :: this

    call print_title('task_manager_type')
    call print("%active: "//this%active)
    call print("%my_worker_id: "//this%my_worker_id)
    call print("%distributed: "//this%distributed)
    call print("%unified_workload: "//this%unified_workload)
    call print("%n_tasks: "//this%n_tasks)
    call print("%n_workers: "//this%n_workers)
    call print_title('MPI_obj')
    call print(this%MPI_obj)
    call print_title('ScaLAPACK_obj')
    call print(this%ScaLAPACK_obj)
    call print_title('task distrubtion')
    call task_manager_show_distribution(this)
    call print_title('')
  end subroutine task_manager_print

end module task_manager_module
