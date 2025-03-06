def run_seeded_simulation(seeded_simulation_parameters_path: str, 
                          experiment_name: str,
                          plot_trajectory: bool) -> None:
    """runs one seeded simulation. This function isolates how to run one simulation 
    from the looping over all seeded simulation parameters of one experiment.
    
    Support for running Slurm requires this function to be importable from its 
    name as reference (simplicity.runners.unit_run.run_seeded_simulation).
    Also this function's arguments are passed by reference (achieved to 
    serializing/deserializing the reference to the simulation parameters 
    that are resolved via the settings_manager.read_seeded_simulation_parameters function).
    """
    # <settings and output managers>
    import simplicity.settings_manager as sm
    import simplicity.output_manager   as om
    parameters       = sm.read_seeded_simulation_parameters(experiment_name, seeded_simulation_parameters_path)
    output_directory = om.setup_output_directory           (experiment_name, seeded_simulation_parameters_path)
    # </settings and output managers>

    ## <simplicity core>        
    import simplicity.simulation       as sim
    simulation = sim.Simplicity       (parameters, output_directory)
    simulation.run()
    ## </simplicity core>        
    if plot_trajectory:    
        simulation.plot()
    
    print(f"run_seeded_simulation: experiment_name={experiment_name}")
