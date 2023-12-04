! Code by: Felipe Reibnitz Willemann
! Computational techniques - Intensive Course Cantabria (2023)
program main
    use tools
    use huckel
    use jacobi
    implicit none
    real(4), allocatable :: coords(:, :), dist(:, :), topol(:, :), H(:, :), U(:, :), pop(:, :), n_occ(:)
    real(4), parameter :: alpha = -11.4, beta = -0.8, error = 1E-10
    real(4) :: total_energy
    character(1), allocatable :: atoms(:)
    character(200) :: file_name
    integer :: i, j, charge, n_electrons

    ! get the name of the input
    write(6, *) "Enter the name of the input file:"
    read(5, '(A)') file_name

    ! get the charge of the molecule
    write(6, *) "Enter the total charge of the molecule"
    read(5, *) charge

    ! read the xyz input file and build the distance matrix
    call read_xyz(file_name, atoms, coords)
    call build_dist(coords, dist)

    ! build the topology matrix and the huckel hamiltonian
    call build_topol(atoms, dist, topol)
    call build_hamiltonian(topol, alpha, beta, H)

    ! diagonalize H (the eigenvectors with be kept in U)
    call diagonalize(H, U, error)

    ! sort H and U
    call sort_eigenvalues(H, U)
    
    ! write the diagonalization results
    write(6, *) "U = "
    call print_matrix(U)
    write(6, *) "H = "
    call print_matrix(H)
    
    ! allocate and initialize the occupancy number vector
    allocate(n_occ(size(H, 1)))
    n_occ = 0.E0
    
    ! define the number of pi electrons
    n_electrons = size(H, 1) - charge

    ! calculate the total energy by adding each electron energy in the MOs
    do i = 1, n_electrons
        ! define the index of the MO to be occupied
        j = 1 + int((i-1)/2)

        ! increase the occupation number
        n_occ(j) = n_occ(j) + 1
        
        ! if the HOMO is degenerated, spins must be parallel (Hund's rule)
        if (abs(H(j, j) - H(j+1, j+1)) < 1.E-4 .and. (n_electrons - i) == 1) then
            n_occ(j+1) = n_occ(j+1) + 1
            exit
        end if        
    end do

    ! print occupation numbers
    write(6, *) "Occupation numbers:"
    write(6, fmt='(*(F12.5))') n_occ

    ! calculate the total energy
    total_energy = sum(matmul(H, n_occ))
    write(6, fmt='(X,A,X,ES12.5,X,A)') "E =", total_energy, "eV"

    ! perform Mulliken analisys
    call mulliken(U, n_occ, pop)
    write(6, *) "Population matrix ="
    call print_matrix(pop)
end program main
