# FitnessCentrality

A non-optimized Julia program to calculate fitness centrality in networks

## Installation

### 1. Install Julia

1. Visit the official Julia website: https://julialang.org/downloads/
2. Download the appropriate version for your operating system
3. Run the installer and follow the prompts

4. Verify the installation:
   - Open a terminal or command prompt
   - Type `julia` and press Enter
   - You should see the Julia REPL (Read-Eval-Print Loop) start up

### 2. Set up the project

1. Open the command line interpreter with:
   ```
   julia --project
   ```
2. Access the package manager by typing:
   ```
   ]
   ```
3. Run the `instantiate` command to download and install the required libraries:
   ```
   instantiate
   ```
4. Exit the package manager by pressing `<backspace>`.
5. Exit `julia` by pressing `<ctrl-d>`.

## Input Files

- Input files contain the list of edges of the network, one per line.
- Nodes can be specified as integers or strings.
- The program automatically detects the number of columns:
  - If more than two columns are found, the network is treated as weighted.
- File extension (e.g., `.net`) is arbitrary; files are essentially CSV or TSV format.

## Running the Program

Execute the program with the following command:

```
julia --project FitnessComplexity.jl filename
```

Where `filename` is the path to the network file you want to analyze. For example:

```
julia --project FitnessComplexity.jl networks/star-wheel.net
```

### Directed Networks

The program will ask if the network is directed. Type `y` if this is the case.

## Results

The analysis results are placed in the `results/` folder.