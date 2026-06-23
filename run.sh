DIR="$1"
LMP="$2"

# Save the current path.
current_dir=$(pwd)

echo "$DIR"
# Move to the target directory to prepare the system there.
cd $DIR
  # Adjust for Your LAMMPS installation and setup of hardware interface.
  case "$LMP" in
    "lmp_mpi")
        mpirun -np 8 lmp_mpi -i run.in.nvt
        ;;
    "lmp_mpi_kokkos")
        lmp_mpi -k on g 1 -sf kk -pk kokkos -in run.in.nvt
        ;;
    "lmp_serial")
        lmp_serial -i run.in.nvt
        ;;
    "lmp_serial_kokkos")
        lmp_serial -k on g 1 -sf kk -pk kokkos -in run.in.nvt
        ;;
    "lmp")
        lmp -i run.in.nvt
        ;;
    "lmp_kokkos")
        lmp -k on g 1 t 20 -sf kk -pk kokkos gpu/aware on neigh half newton off -in run.in.nvt
        ;;
    *)
        echo "Unknown LAMMPS run setup."
        ;;
  esac

# Return to the initial user path.
cd current_dir
