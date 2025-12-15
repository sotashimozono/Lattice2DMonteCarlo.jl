# --- Rules ---
const metropolis = Metropolis()
const glauber = Glauber()
export metropolis, glauber

# --- Selections ---
const random_select = RandomSiteSelection()
const sequential_sweep = SequentialSweep()
export random_select, sequential_sweep

# --- Proposals ---
const spin_flip = SpinFlip()
const spin_exchange = SpinExchange()
const uniform_shift = UniformShift()
export spin_flip, spin_exchange, uniform_shift

# --- Algorithms ---
const localupdate = LocalUpdate()
export localupdate
