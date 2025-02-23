using Random, Plots
using FFTW

# Define parameters
const L = 50  # Grid size
const P_init = 0.3  # Initial prey density
const Q_init = 0.2  # Initial predator density
r_P = 0.7  # Prey reproduction probability
r_Q = 0.7  # Predators preying sucsyfully probability
const T_starve = 1  # Predator starvation time
steps = 1000  # Number of simulation steps

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
            if grid[ni, nj] == PREY  # Eat prey
                new_grid[ni, nj] = PREDATOR
                new_starvation[ni, nj] = 0  # Reset starvation at new position
            elseif new_grid[ni, nj] == EMPTY
                new_grid[ni, nj] = PREDATOR
                new_starvation[ni, nj] = starvation_time[i, j] + 1
                new_starvation[i, j] = 0  # Move predator
                new_grid[i, j] = EMPTY  # Clear old predator position
            end

            if new_starvation[ni, nj] > T_starve 
                new_grid[ni, nj] = EMPTY  # Predator dies at its new location due to starvation
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

# Run simulation
grid, starvation_time = initialize_grid(L, P_init, Q_init)
prey_counts = Int[]
predator_counts = Int[] 

anim = @animate for step in 1:steps
    update_grid!(grid, starvation_time, L, r_P, T_starve)
    heatmap(grid, c=:viridis, title="Prey-Predator Model on $(L)×$(L) Lattice- Step $step", size=(800,750))
    
    # Count and store populations
    num_prey, num_predators = count_species(grid)
    push!(prey_counts, num_prey)
    push!(predator_counts, num_predators)
end

# Save animation
gif(anim, "prey_predator_simulation.gif", fps=100)

# Plot prey and predator populations over time
 plot(1:steps, [prey_counts, predator_counts], 
    label=["Prey" "Predators"], 
    xlabel="Time Step", ylabel="Population",
    title="Prey-Predator Population Dynamics",
    lw=2, legend=:topright, size=(800,600))