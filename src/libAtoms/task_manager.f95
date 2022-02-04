module task_manager_module
  !> A task manager determines the distribution of tasks to workers with
  !> respect to the tasks' demands.
  !>
  !> Each task stores integer data in %idata.
  !> %idata(1) is used for distribution.
  !>
  !> The special %worker_id=SHARED is used for shared tasks. Each worker
  !> has a local entry with local %idata. %idata(1) contains the local
  !> size, %idata(2) contains the global offset.

  use libatoms_module, only : initialise, finalise, print_title, print, &
      operator(//), system_abort, inoutput, optional_default, OUTPUT, PRINT_VERBOSE
  use linearalgebra_module, only : heap_sort
  use MPI_context_module, only : MPI_Context, print
  use ScaLAPACK_module, only : Scalapack, print

  implicit none
  private

  public :: task_manager_type, task_type, worker_type, SHARED
  public :: task_manager_init_tasks, task_manager_init_workers, task_manager_deallocate
  public :: task_manager_init_idata
  public :: task_manager_add_task, task_manager_add_worker
  public :: task_manager_distribute_tasks, task_manager_check_distribution
  public :: task_manager_print, task_manager_export_distribution

  integer, parameter :: UNDEFINED = -1
  integer, parameter :: SHARED = 0

  type task_type
    integer :: index = UNDEFINED !> initial task index
    integer, dimension(:), allocatable :: idata !> integer data for this task
    integer :: worker_id = UNDEFINED !> ID of worker assigned to this task
  end type task_type

  type worker_type
    integer :: index = UNDEFINED !> initial worker index
    integer :: n_tasks = 0 !> number of tasks assigned to this worker
    integer :: n_padding = 0 !> number of padding rows to reach unified workload
    type(task_type), dimension(:), allocatable :: tasks !> tasks assigned to this worker
  end type worker_type

  type task_manager_type
    logical :: active = .false. !> whether task manager is used
    logical :: distributed = .false. !> whether tasks have been distributed
    integer :: my_worker_id = 0 !> compute only tasks for this worker
    integer :: unified_workload = 0 !> highest workload among workers, target for padding

    integer :: n_tasks = 0 !> number of tasks
    integer :: n_workers = 0 !> number of workers
    integer :: n_shared = 0 !> number of shared tasks

    integer, dimension(:), allocatable :: idata !< custom shared integer data

    type(MPI_Context) :: MPI_obj
    type(ScaLAPACK) :: ScaLAPACK_obj

    type(task_type), dimension(:), allocatable :: tasks !> tasks to compute
    type(worker_type), dimension(:), allocatable :: workers !> workers to compute tasks
  end type task_manager_type

  contains

  !> Add task to list
  subroutine task_manager_add_task(this, idata1, n_idata, worker_id)
    type(task_manager_type), intent(inout) :: this
    integer, intent(in), optional :: idata1 !> first idata entry for new task
    integer, intent(in), optional :: n_idata !> number of idata entries for new task
    integer, intent(in), optional :: worker_id !> preset worker for this task

    integer :: my_n_idata

    if (.not. this%active) return
    my_n_idata = optional_default(1, n_idata)

    this%n_tasks = this%n_tasks + 1
    if (this%n_tasks > size(this%tasks)) call system_abort("More tasks added than allocated.")
    associate(task => this%tasks(this%n_tasks))
        allocate(task%idata(my_n_idata))
        task%idata = UNDEFINED
        task%index = this%n_tasks
        if (present(idata1)) task%idata(1) = idata1
        if (present(worker_id)) task%worker_id = worker_id
    end associate
  end subroutine task_manager_add_task

  !> Add worker to list
  subroutine task_manager_add_worker(this)
    type(task_manager_type), intent(inout) :: this

    if (.not. this%active) return

    this%n_workers = this%n_workers + 1
    if (this%n_workers > size(this%workers)) call system_abort("More workers added than allocated.")
    associate(worker => this%workers(this%n_workers))
        worker%index = this%n_workers
        worker%n_tasks = 0
    end associate
  end subroutine task_manager_add_worker

  !> @brief Distribute tasks among workers for minimal idata(1) usage
  subroutine task_manager_distribute_tasks(this)
    type(task_manager_type), intent(inout) :: this

    integer :: i, t, w ! indices
    integer :: n
    integer, dimension(:), allocatable :: workloads
    integer, dimension(:), allocatable :: idata1_list
    integer, dimension(:), allocatable :: index_list

    if (.not. this%active) return
    if (this%n_tasks < 1) return
    if (this%n_workers < 1) call system_abort("task_manager_distribute_tasks: No workers initialised.")

    allocate(workloads(this%n_workers))
    allocate(idata1_list(this%n_tasks))
    allocate(index_list(this%n_tasks))

    ! copy to temp arrays, sort both by idata1
    do i = 1, this%n_tasks
      if (.not. allocated(this%tasks(i)%idata)) call system_abort("task_manager_distribute_tasks: idata of task "//i//" is unallocated.")
      idata1_list(i) = this%tasks(i)%idata(1)
      index_list(i) = this%tasks(i)%index
    end do
    call heap_sort(idata1_list, i_data=index_list) ! index_list sorted via idata1_list

    ! count workload of preset tasks
    workloads = 0
    do t = 1, this%n_tasks
        i = index_list(t)
        w = this%tasks(i)%worker_id
        if (w >= 1) then
          this%workers(w)%n_tasks = this%workers(w)%n_tasks + 1
          workloads(w) = workloads(w) + idata1_list(t)
        end if
    end do

    ! distribute tasks to workers to minimize maximum idata1 sum
    do t = this%n_tasks, 1, -1
        i = index_list(t)
        if (this%tasks(i)%worker_id /= UNDEFINED) cycle ! skip preset tasks
        w = minloc(workloads, 1)
        this%tasks(i)%worker_id = w
        this%workers(w)%n_tasks = this%workers(w)%n_tasks + 1
        workloads(w) = workloads(w) + idata1_list(t)
    end do

    this%unified_workload = maxval(workloads)

    ! fill workers with task id lists (replication for simpler retrieval)
    this%n_shared = count(this%tasks(:)%worker_id == SHARED)
    do w = 1, this%n_workers
        n = this%workers(w)%n_tasks + this%n_shared
        allocate(this%workers(w)%tasks(n))
        this%workers(w)%n_tasks = 0
    end do
    do t = 1, this%n_tasks
        w = this%tasks(t)%worker_id
        if (w < 1) cycle
        n = this%workers(w)%n_tasks + 1
        this%workers(w)%tasks(n) = this%tasks(t)
        this%workers(w)%n_tasks = n
    end do

    do w = 1, this%n_workers
      this%workers(w)%n_padding = this%unified_workload - workloads(w)
    end do

    call task_manager_distribute_shared_tasks(this)

    this%distributed = .true.
  end subroutine task_manager_distribute_tasks

  subroutine task_manager_distribute_shared_tasks(this)
    ! fill padding first, then distribute rows evenly
    type(task_manager_type), intent(inout) :: this
  
    integer :: i, n, t, w
    integer :: rest, n_split, n_extra, n_plus, offset

    if (.not. this%active) return

    do t = 1, this%n_tasks
      if (this%tasks(t)%worker_id /= SHARED) cycle

      rest = max(this%tasks(t)%idata(1) - sum(this%workers(:)%n_padding), 0)
      n_split = rest / this%n_workers
      n_extra = mod(rest, this%n_workers)
      n_plus = n_split
      if (n_extra > 0) n_plus = n_split + 1

      rest = this%tasks(t)%idata(1)
      offset = 0
      do w = 1, this%n_workers
        n = this%workers(w)%n_tasks + 1
        associate(worker => this%workers(w), w_task => this%workers(w)%tasks(n))
          call task_deepcopy(w_task, this%tasks(t))

          w_task%idata(1) = worker%n_padding + n_split
          if (w <= n_extra) w_task%idata(1) =  w_task%idata(1) + 1
          w_task%idata(1) = min(w_task%idata(1), rest)
          w_task%idata(2) = offset
          worker%n_padding = worker%n_padding + n_plus - w_task%idata(1)

          rest = rest - w_task%idata(1)
          offset =  offset + w_task%idata(1)
          worker%n_tasks = n
        end associate
      end do
      this%unified_workload = this%unified_workload + n_plus
    end do
  end subroutine task_manager_distribute_shared_tasks

  subroutine task_manager_check_distribution(this)
    type(task_manager_type), intent(in) :: this

    integer :: t, w, n, tmp

    if (.not. this%active) return
    if (.not. this%distributed) return

    tmp = task_manager_get_max_idata1_sum(this)
    if (tmp /= this%unified_workload) call system_abort("task_manager_check_distribution: unified_workload is inconsistent: "//tmp//" "//this%unified_workload)
    if (any(this%workers(:)%n_padding < 0)) call system_abort("task_manager_check_distribution: negative padding: "//this%workers(:)%n_padding)
  end subroutine task_manager_check_distribution

  subroutine task_manager_show_distribution(this)
    type(task_manager_type), intent(in) :: this

    integer :: i, t, w, idata1_sum, n_tasks

    if (.not. this%active) return

    do w = 1, this%n_workers
        call print("Tasks of worker "//w)
        idata1_sum = 0
        n_tasks = 0
        do i = 1, this%workers(w)%n_tasks
          t = this%workers(w)%tasks(i)%index
          call print(i // ": " // t // ": " // this%workers(w)%tasks(i)%idata)
          idata1_sum = idata1_sum + this%tasks(t)%idata(1)
        end do
        call print("Summary: "//this%workers(w)%n_tasks//" tasks, "//idata1_sum//" idata1_sum, "//this%workers(w)%n_padding//" padding")
    end do
  end subroutine task_manager_show_distribution

  subroutine task_manager_export_distribution(this, only_id)
    type(task_manager_type), intent(in) :: this
    integer, optional :: only_id !> only this worker writes distribution

    integer :: i, t, w
    integer :: only_id_opt
    type(Inoutput) :: file

    if (.not. this%active) return
    only_id_opt = optional_default(UNDEFINED, only_id)
    if (only_id_opt >= 0 .and. only_id_opt /= this%my_worker_id) return

    call initialise(file, "tm.dist", OUTPUT)
    call print("worker_id task_id_local task_id_global idata", file=file)
    do w = 1, this%n_workers
        do i = 1, this%workers(w)%n_tasks
          t = this%workers(w)%tasks(i)%index
          call print(w // " " // i // " " // t // " " // this%workers(w)%tasks(i)%idata, file=file)
        end do
    end do
    call finalise(file)
  end subroutine task_manager_export_distribution

  subroutine task_manager_init_tasks(this, n_tasks)
    type(task_manager_type), intent(inout) :: this
    integer, intent(in) :: n_tasks !> number of tasks to allocate

    if (.not. this%active) return

    call print("Allocate task_manager%tasks for "//n_tasks//" tasks.", PRINT_VERBOSE)
    allocate(this%tasks(n_tasks))
  end subroutine task_manager_init_tasks

  subroutine task_manager_init_workers(this, n_workers)
    type(task_manager_type), intent(inout) :: this
    integer, intent(in) :: n_workers !> number of workers to allocate and set

    integer :: w

    if (.not. this%active) return

    allocate(this%workers(n_workers))
    do w = 1, n_workers
      call task_manager_add_worker(this)
    end do
  end subroutine task_manager_init_workers

  subroutine task_manager_init_idata(this, n_idata)
    type(task_manager_type), intent(inout) :: this
    integer, intent(in) :: n_idata

    if (.not. this%active) return

    allocate(this%idata(n_idata))
  end subroutine task_manager_init_idata

  subroutine task_manager_deallocate(this)
    type(task_manager_type), intent(inout) :: this

    this%distributed = .false.
    this%n_tasks = 0
    this%n_workers = 0
    if (allocated(this%idata)) deallocate(this%idata)
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

  function task_manager_get_max_idata1_sum(this) result(res)
    type(task_manager_type), intent(in) :: this

    integer :: res
    integer :: n, t, w

    if (.not. this%active) return

    res = 0
    do w = 1, this%n_workers
      n = sum([(this%workers(w)%tasks(t)%idata(1), t = 1, this%workers(w)%n_tasks)])
      res = max(n, res)
    end do
  end function task_manager_get_max_idata1_sum

  subroutine task_deepcopy(this, from)
    type(task_type), intent(inout) :: this
    type(task_type), intent(in) :: from

    this%index = from%index
    allocate(this%idata, source=from%idata)
    this%worker_id = from%worker_id
  end subroutine task_deepcopy

end module task_manager_module
