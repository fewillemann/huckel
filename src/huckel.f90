! Code by: Felipe Reibnitz Willemann
! Computational techniques - Intensive Course Cantabria (2023)
module huckel
    contains
    subroutine read_xyz(file_name, atoms, coords)
        implicit none
        character(100), intent(in) :: file_name
        character(len=1), allocatable, intent(out) :: atoms(:)
        real(4), allocatable, intent(out) :: coords(:, :)
        integer :: i, n_atoms, un
        
        ! first we open the xyz input to read the first line (number of atoms) and skip the second (comment)
        open(newunit=un, file=trim(file_name), status="old")
        read(un, *) n_atoms
        read(un, *) 
        
        ! allocate the arrays acording to the number os atoms
        allocate(atoms(n_atoms), coords(n_atoms, 3))

        ! run through the rest of the file
        do i = 1, n_atoms
            ! assign each atoms name and coordinates to the arrays
            read(un, *) atoms(i), coords(i, :)
        end do
        close(un)

        return
    end subroutine read_xyz
    
    subroutine build_dist(coords, dist)
        implicit none
        real(4), intent(in) :: coords(:, :)
        real(4), allocatable, intent(out) :: dist(:, :)
        real(4) :: diff(3)
        integer :: i, j, n_atoms

        ! get the number of atoms
        n_atoms = size(coords, 1)
        
        ! allocate and initialize the distante matrix
        allocate(dist(n_atoms, n_atoms))
        dist = 0.E0

        ! run through the lower half of the distante matrix (diagonal is zero)
        do i = 2, n_atoms
            do j = 1, i-1
                ! calculate the difference vector between atoms i and j
                diff = coords(i, :) - coords(j, :)

                ! calculate the distance (norm of the difference)
                dist(i, j) = norm2(diff)

                ! the distance matrix is symmetric
                dist(j, i) = dist(i, j)
            end do
        end do

        return
    end subroutine build_dist

    subroutine build_topol(atoms, dist, topol)
        implicit none
        character(1), intent(in) :: atoms(:)
        real(4), intent(in) :: dist(:, :)
        real(4), allocatable, intent(out) :: topol(:, :)
        real(4), parameter :: dist_C = 1.54
        integer, allocatable :: idx_C(:)
        integer :: i, j, n_atoms, n_C

        ! get the number of atoms
        n_atoms = size(atoms, 1)

        ! initiate the number of carbons and allocate the carbons indexes array (in the distance matrix)
        n_C = 0
        allocate(idx_C(n_atoms))
        
        ! get number and indexes of the carbons
        do i = 1, n_atoms
            if (atoms(i) == 'C') then
                n_C = n_C + 1
                idx_C(n_C) = i
            end if
        end do

        ! allocate and initialize the topology matrix
        allocate(topol(n_C, n_C))
        topol = 0.E0

        ! run through the lower half of the topol matrix
        do i = 2, n_C
            do j = 1, i-1
                ! check if the distance between the ith and jth carbons allows a bond
                if (dist(idx_C(i), idx_C(j)) < dist_C) then
                    ! carbons are connected 
                    topol(i, j) = 1.E0

                    ! topology matrix is symmetric
                    topol(j, i) = 1.E0
                end if
            end do
        end do

        deallocate(idx_C)

        return
    end subroutine build_topol

    subroutine build_hamiltonian(topol, alpha, beta, H)
        implicit none
        real(4), intent(in) :: topol(:, :), alpha, beta
        real(4), allocatable, intent(out) :: H(:, :)
        integer :: i, n_C

        ! get the number of carbons (dimension of the huckel hamiltonian)
        n_C = size(topol, 1)

        ! allocate and initialize the hamiltonian as the beta matrix (conected carbons)
        allocate(H(n_C, n_C))
        H = topol * beta

        ! define the diagonal elements (alpha)
        do i = 1, n_C
            H(i, i) = alpha
        end do
        
        return
    end subroutine build_hamiltonian

    subroutine mulliken(U, n_occ, pop)
        implicit none
        real(4), intent(in) :: U(:, :), n_occ(:)
        real(4), allocatable, intent(out) :: pop(:, :)
        integer :: i, j, n_atoms
    
        ! get the number of atoms (pz atomic orbitals)
        n_atoms = size(U, 1)
    
        ! allocate and initialize the population matrix
        allocate(pop(n_atoms, n_atoms))
        pop = 0.E0
    
        ! run through the ith atomic orbitals
        do i = 1, n_atoms
            ! run through the jth atomic orbitals
            do j = 1, i
                ! perform the Mulliken sum 
                pop(i, j) = sum(n_occ(:) * U(i, :) * U(j, :))

                ! the matrix is symmetrix
                pop(j, i) = pop(i, j)
            end do
        end do
    
        return
    end subroutine mulliken
end module huckel
