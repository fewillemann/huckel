! Code by: Felipe Reibnitz Willemann
! Computational techniques - Intensive Course Cantabria (2023)
module jacobi
    contains
    subroutine rotate(s, c, p, q, A, U)
        implicit none
        real(4), intent(in) :: c, s
        real(4), intent(out) :: A(:, :), U(:, :)
        real(4) :: elp, elq, elr
        integer, intent(in) :: p, q
        integer :: r

        ! keep the old diagonal values
        elp = A(p, p)
        elq = A(q, q)

        ! transform the diagonal of A
        A(p, p) = (c**2) * elp - (2*s*c) * A(p, q) + (s**2) * elq
        A(q, q) = (s**2) * elp + (2*s*c) * A(p, q) + (c**2) * elq

        ! run through the lines of the matrices
        do r = 1, size(A, 1)
            ! diagonal elements of A already transformed
            if((r /= p) .and. (r /= q)) then
                ! keep the old value
                elr = A(r, p)

                ! calculate new off-diagonal elements
                A(r, p) = c * elr - s * A(r, q)
                A(r, q) = s * elr + c * A(r, q)
                
                ! the matrix is symmetric
                A(p, r) = A(r, p)
                A(q, r) = A(r, q)
            end if

            ! keep the old value
            elr = U(r, p)

            ! rotate U
            U(r, p) = c * elr - s * U(r, q)
            U(r, q) = s * elr + c * U(r, q)
        end do

        return
    end subroutine rotate

    subroutine calculate_delta(A, delta)
        implicit none
        real(4), intent(in) :: A(:, :)
        real(4), intent(out) :: delta
        integer :: i, j

        ! initialize delta
        delta = 0

        ! run through the lower half of A
        do i = 2, size(A, 1)
            do j = 1, i-1
                ! calculate delta by summing up the squares
                delta = delta + A(i, j) ** 2
            end do
        end do

        ! take the square root
        delta = sqrt(delta)

        return
    end subroutine calculate_delta

    subroutine diagonalize(A, U, error)
        implicit none
        real(4), intent(in) :: error
        real(4), intent(out) :: A(:, :)
        real(4), allocatable, intent(out) :: U(:, :)
        real(4) :: s, c, t, sigma, delta
        integer :: p, q, n
        
        ! get A's dimension
        n = size(A, 1)

        ! allocate and initialize U (eigenvectors matrix) as the identity
        allocate(U(n, n))
        U = 0.E0
        do p = 1, n
            U(p, p) = 1.E0
        end do

        ! initialize delta and error
        delta = 1.E0

        ! Repeat swipes until convergence
        do while (delta > error)
            ! iterate over the lower half of the matrix (swipe)
            do p = 2, n
                do q = 1, p-1
                    ! skip the calculation if the target element is already zero
                    if (A(p, q) .eq. 0.E0) cycle

                    ! otherwise, calculate the rotation parameter
                    sigma = (A(q, q) - A(p, p)) / (2.E0 * A(p, q))

                    ! calculate the smallest value of the tangent
                    t = sign(1.E0, sigma) / (abs(sigma) + sqrt(sigma ** 2 + 1))

                    ! and the other trig functions
                    c = 1.E0 / (sqrt(1.E0 + t**2))
                    s = t * c

                    ! apply the (p, q) givens rotation to A and U
                    call rotate(s, c, p, q, A, U)

                    ! zero the target elements
                    A(p, q) = 0.E0
                    A(q, p) = 0.E0
                end do
            end do
            ! after the swipe, calculate delta to check convergence
            call calculate_delta(A, delta)
        end do

        ! retrieve just the diagonal
        do p = 2, n
            do q = 1, p-1
                A(p, q) = 0.E0
                A(q, p) = 0.E0
            end do
        end do
        
        return
    end subroutine diagonalize
end module jacobi