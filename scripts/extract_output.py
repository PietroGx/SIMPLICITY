import os
import shutil
import pandas as pd
import random 
import simplicity.config as config 

def process_experiment_folder(experiment_name):
    """
    Process the main experiment folder to extract .png and .csv files from its immediate subdirectories.
    The files are copied and renamed before being stored in the '05_summary' directory.

    Args:
        experiment_name (str): The name of the experiment whose files are to be processed.
    """
    # Set up paths
    data_dir = config.get_data_dir()
    experiment_dir = os.path.join(data_dir, experiment_name)
    summary_dir = os.path.join(experiment_dir, '05_summary')

    # Create the summary directory if it doesn't exist
    if not os.path.exists(summary_dir):
        os.makedirs(summary_dir)

    # Iterate over each subdirectory in the experiment folder
    for subdir in os.listdir(experiment_dir):
        subdir_path = os.path.join(experiment_dir, subdir)
        if os.path.isdir(subdir_path):
            # Process files directly in the subdirectory
            for file in os.listdir(subdir_path):
                if file.endswith('.png') or file.endswith('.csv'):
                    file_path = os.path.join(subdir_path, file)
                    new_filename = f'{os.path.splitext(file)[0]}_{subdir}{os.path.splitext(file)[1]}'
                    shutil.copy(file_path, os.path.join(summary_dir, new_filename))

def process_trajectory_files(experiment_name):
    """
    Process seed subdirectories within experiment_folder/04_Output/simulations 
    to select the 3 longest simulations.
    Extract the 'simulation_trajectory.png' from these selected seeds and store 
    them in the '05_summary' directory.

    Args:
        experiment_name (str): The name of the experiment whose trajectory files are to be processed.
    """
    # Set up paths
    data_dir = config.get_data_dir()
    experiment_dir = os.path.join(data_dir, experiment_name)
    output_dir = os.path.join(experiment_dir,'04_Output')
    summary_dir = os.path.join(experiment_dir, '05_summary')
    
    # Iterate over each subfolder in the experiment output folder
    for subdir in os.listdir(output_dir):
        subdir_path = os.path.join(output_dir, subdir)
        
        # Gather seed_number directories and their trajectory.csv last time entries
        seed_folders = []
        if os.path.isdir(subdir_path):
            # Iterate over each seed folder in the subfolder
            for seed_dir in os.listdir(subdir_path):
                seed_path = os.path.join(subdir_path, seed_dir)
                
                trajectory_file = os.path.join(seed_path, 'trajectory.csv')
                if os.path.exists(trajectory_file):
                    df = pd.read_csv(trajectory_file)
                    last_time = df['time'].iloc[-1]
                    seed_folders.append((last_time, seed_path, seed_dir))
                    
            # Sort seed folders by last time entry and select top 3, random in case of ties
            seed_folders.sort(reverse=True, key=lambda x: x[0])
            top_3_seed_folders = sorted(seed_folders[:3], key=lambda x: random.random())
    
            # Create summary subfolder directory 
            save_subdir = os.path.join(summary_dir, subdir)
            if not os.path.exists(save_subdir):
                os.makedirs(save_subdir)

        # Copy and rename simulation_trajectory.png for each selected seed folder
        for _, seed_path, seed_dir in top_3_seed_folders:
            sim_traj_src = os.path.join(seed_path, 'simulation_trajectory.png')
            if os.path.exists(sim_traj_src):
                sim_traj_dest = os.path.join(save_subdir, f'simulation_trajectory_{seed_dir}.png')
                shutil.copy(sim_traj_src, sim_traj_dest)


if __name__ == "__main__":
    # if len(sys.argv) != 2:
    #     print("Usage: python script.py <experiment_name>")
    #     sys.exit(1)

    # experiment_name = sys.argv[1]
    experiment_name = 'test_refactoring'
    process_experiment_folder(experiment_name)
    process_trajectory_files(experiment_name)

