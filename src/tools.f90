! Code by: Felipe Reibnitz Willemann
! Computational techniques - Intensive Course Cantabria (2023)
module tools
    implicit none
    contains
    subroutine print_matrix(mat)
        implicit none      
        real(4), intent(in) :: mat(:, :)
        integer :: i
		
		! run through the matrix, writing every element of each line per instance
        do i = 1, size(mat, 1)
            write(6, fmt='(*(F12.5))') mat(i, :)
        end do

        return
    end subroutine print_matrix

    subroutine sort_eigenvalues(D, U)
        implicit none
        real(4), intent(out) :: D(:, :), U(:, :)
        real(4), allocatable :: gr_vec(:)
        real(4) :: gr_val
        integer :: i, j, n 

		! get the dimension of the matrices
        n = size(D, 1)

		! allocate the vector to keep the greater eigenvector
        allocate(gr_vec(n))

		! start the bubble loop
        do i = 1, n-1
            do j = 1, n-1
                if (D(j, j) > D(j+1, j+1)) then
			        ! keep the greater eigenvalue and eigenvector
                    gr_val = D(j+1, j+1)
                    gr_vec = U(:, j+1)

			        ! switch eigenvalue positions
                    D(j+1, j+1) = D(j, j)
                    D(j, j) = gr_val

			        ! switch eigenvectors positions
                    U(:, j+1) = U(:, j)
                    U(:, j) = gr_vec
                end if
            end do
        end do

        return
    end subroutine sort_eigenvalues
end module tools
