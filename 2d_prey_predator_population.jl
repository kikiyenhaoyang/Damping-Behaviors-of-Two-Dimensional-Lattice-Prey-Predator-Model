using Random, Plots
using Statistics

# Define parameters
L = 50  # Grid size
P_init = 0.3  # Initial prey density
Q_init = 0.2  # Initial predator density
r_P = 0.9  # Prey reproduction probability
r_Q = 0.7  # Predators preying sucsyfully probability
T_starve = 1  # Predator starvation time
steps = 1000  # Number of simulation steps
simulations = 5000  # Number of simulations

# Define states
const EMPTY = 0
const PREY = 1
const PREDATOR = 2

# Initialize grid
function initialize_grid(L, P_init, Q_init)
    grid = fill(EMPTY, L, L)
    starvation_time = fill(0, L, L)  # Track starvation time for predators

    for i in 1:L, j in 1:L
        r = rand()
        if r < P_init
            grid[i, j] = PREY
        elseif r < P_init + Q_init
            grid[i, j] = PREDATOR
        end
    end
    return grid, starvation_time
end

# Get random neighborhood moves
const MOVES = [(0,1), (1,0), (0,-1), (-1,0), (-1,-1), (-1,1), (1,-1), (1,1)]

function get_neighbors(i, j, L)
    random_move = MOVES[rand(1:8)]
    ni, nj = i + random_move[1], j + random_move[2]  # Random move
    return (mod1(ni, L), mod1(nj, L))  # Periodic boundary
end

# Update grid for one Monte Carlo step
function update_grid!(grid, starvation_time, L, r_P, T_starve)
    new_grid = copy(grid)
    new_starvation = copy(starvation_time)

    for _ in 1:L*L  # Process L^2 random sites
        i, j = rand(1:L), rand(1:L)
        ni, nj = get_neighbors(i, j, L)

        if grid[i, j] == PREY
            if new_grid[ni, nj] == EMPTY && rand() < r_P
                new_grid[ni, nj] = PREY  # Reproduce
            end

        elseif grid[i, j] == PREDATOR
            if grid[ni, nj] == PREY && rand() < r_Q # Eat prey
                new_grid[ni, nj] = PREDATOR
                new_starvation[ni, nj] = 0  # Reset starvation at new position
            elseif new_grid[ni, nj] == EMPTY
                new_grid[ni, nj] = PREDATOR
                new_starvation[ni, nj] = starvation_time[i, j] + 1
                new_starvation[i, j] = 0  # Move predator
                new_grid[i, j] = EMPTY  # Clear old predator position
            end

            if new_starvation[ni, nj] > T_starve
                new_grid[ni, nj] = EMPTY  # Predator dies at its new location
            end
        end
    end

    grid .= new_grid
    starvation_time .= new_starvation
end

# Count the number of prey and predators in the grid
function count_species(grid)
    num_prey = count(x -> x == PREY, grid)
    num_predators = count(x -> x == PREDATOR, grid)
    return num_prey, num_predators
end

# Run multiple simulations and compute averages
prey_totals = zeros(steps)  # Array to store total prey counts
predator_totals = zeros(steps)  # Array to store total predator counts

for t in 1:simulations
    grid, starvation_time = initialize_grid(L, P_init, Q_init)
    
    for step in 1:steps
        update_grid!(grid, starvation_time, L, r_P, T_starve)
        
        # Count populations
        num_prey, num_predators = count_species(grid)
        prey_totals[step] += num_prey
        predator_totals[step] += num_predators
    end
end

# Compute the average prey and predator populations per step
avg_prey = prey_totals ./ simulations
avg_predators = predator_totals ./ simulations

# Plot the averaged population dynamics
plot(1:steps, [avg_prey / ( L * L ), avg_predators / ( L * L )], 
    label=["Avg Prey" "Avg Predators"], 
    xlabel="Time Step", ylabel="Average Population",
    title="Average Prey-Predator Population Density
    of $(L)×$L Lattices Over $simulations Simulations",
    lw=2, legend=:topright)

# Compute the FFT
fft_signal_prey, fft_signal_predator = fft(avg_prey / ( L * L ) .- mean(avg_prey[501:1000])/(L*L)), fft(avg_predators / ( L * L ) .- mean(avg_predators[501:1000]) / ( L * L ))

# Compute the frequencies
N = steps
freqs = fftfreq(N, 1/N)  # Frequency bins

# Plot the magnitude spectrum
plot(freqs[1:N÷2], abs.(fft_signal_prey[1:N÷2]),
    xlabel="Frequency (Hz)", ylabel="Amplitude",
    title="Fourier Transform of Prey Population Density
    on $(L)×$L Lattices Over $simulations Simulations",
    lw=2, legend=false)

plot(freqs[1:N÷2], abs.(fft_signal_predator[1:N÷2]),
    xlabel="Time Step", ylabel="Average Population",
    title="Fourier Transform of Predator Population Density
    on $(L)×$L Lattices Over $simulations Simulations",
    lw=2, legend=false)