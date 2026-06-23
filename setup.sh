DIR="$1"

# Save the current path.
current_dir=$(pwd)

# Move to the target directory to prepare the system there.
cd $DIR
  # Assemble the system with defined cell size and initial positions of
  # polymers.
  moltemplate.sh system.lt

  # Remove the temporal and guide files that will not be used for simulations.
  rm -rf output_ttree/
  rm run.in.EXAMPLE

# Return to the initial user path.
cd $current_dir
