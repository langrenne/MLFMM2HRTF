! Copyright (c) 2024-2026 Conservatoire national des arts et metiers (Cnam), Christophe Langrenne
! SPDX-License-Identifier: MIT
module unique_sort_mod
! Porte depuis Tuyaux/.../modules_sort.f08 (W. M. Peterson, domaine public).
! Tri + unicite d'un tableau d'entiers. Utilise par lecture_partitions_ver4 (T3).
    contains

    subroutine unique_sort_allocatable_arg(val, n, final)
        implicit none
        integer :: i, n
        integer :: min_val, max_val
        integer, dimension(n) :: val, unique
        integer, dimension(:), allocatable, intent(inout) :: final

        i = 0
        min_val = minval(val) - 1
        max_val = maxval(val)
        do while (min_val < max_val)
            i = i + 1
            min_val = minval(val, mask=val > min_val)
            unique(i) = min_val
        enddo

        if (allocated(final)) then
            final = unique(1:i)
        else
            allocate(final(i), source=unique(1:i))
        endif
    end subroutine unique_sort_allocatable_arg

    subroutine unique_sort(val, n, unique, i)
        implicit none
        integer :: i, n
        integer :: min_val, max_val
        integer, dimension(n) :: val, unique

        i = 0
        min_val = minval(val) - 1
        max_val = maxval(val)
        do while (min_val < max_val)
            i = i + 1
            min_val = minval(val, mask=val > min_val)
            unique(i) = min_val
        enddo
    end subroutine unique_sort

end module unique_sort_mod
