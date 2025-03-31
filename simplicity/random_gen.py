def randomgen(seed=None):
    import numpy as np
    # return instance of random number generator
    if seed == None:
        raise Exception("DEPRECATED: seed must be set.")
        return np.random.default_rng()        
    else:
        return np.random.default_rng(seed)

    
