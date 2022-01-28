program test_task_manager
  use system_module
  use task_manager_module

  call system_initialise()
  call test_1to9_for_2()
  call print("test_task_manager OK")
  call system_finalise()

  contains

  subroutine test_1to9_for_2()
    integer, parameter :: n_workers = 2
    integer, parameter :: idata1_normal(*) = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    integer, parameter :: idata1_shared(*) = [10]

    integer, parameter :: workers_n_tasks(*) = [6, 5]
    integer, parameter :: worker1_task_ids(*) = [1, 2, 5, 6, 9, 10]
    integer, parameter :: worker2_task_ids(*) = [3, 4, 7, 8, 10]
    integer, parameter :: workers_paddings(*) = [0, 1]
    integer, parameter :: unified_workload = 28

    type(task_manager_type) :: task_manager

    call test_core(task_manager, n_workers, idata1_normal, idata1_shared)
    if (any(task_manager%workers(:)%n_tasks /= workers_n_tasks)) call system_abort("test_1to9_for_2:workers_n_tasks FAIL")
    if (any(task_manager%workers(1)%tasks(:)%index /= worker1_task_ids)) call system_abort("test_1to9_for_2:worker1_task_ids FAIL")
    if (any(task_manager%workers(2)%tasks(:)%index /= worker2_task_ids)) call system_abort("test_1to9_for_2:worker2_task_ids FAIL")
    if (task_manager%unified_workload /= unified_workload) call system_abort("test_1to9_for_2:unified_workload FAIL")
    if (any(task_manager%workers(:)%n_padding /= workers_paddings)) call system_abort("test_1to9_for_2:workers_paddings FAIL")
  end subroutine test_1to9_for_2

  subroutine test_core(task_manager, n_workers, idata1_normal, idata1_shared)
    integer, intent(in) :: n_workers
    integer, intent(in) :: idata1_normal(:), idata1_shared(:)
    
    type(task_manager_type), intent(inout) :: task_manager

    integer :: t

    task_manager%active = .true.
    call task_manager_init_workers(task_manager, n_workers)
    call task_manager_init_tasks(task_manager, size(idata1_normal) + size(idata1_shared))

    do t = 1, size(idata1_normal)
      call task_manager_add_task(task_manager, idata1_normal(t))
    end do

    do t = 1, size(idata1_shared)
      call task_manager_add_task(task_manager, idata1_shared(t), n_idata=2, worker_id=SHARED)
    end do

    call task_manager_distribute_tasks(task_manager)
  end subroutine test_core

end program test_task_manager
