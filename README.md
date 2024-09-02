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
5. Exit `julia` by pressing `<ctrl+d>`.


### 2. Set up the project

After installing Julia, you need to set up the project environment. This process ensures that all required dependencies are installed and your project is isolated from other Julia projects.

1. Open a terminal or command prompt.

2. Navigate to the directory containing the FitnessCentrality project files.

3. Start Julia in project mode by typing:
   ```
   julia --project
   ```
   This command starts Julia and activates the project environment defined in the current directory.

4. Once Julia starts, you'll see the Julia REPL (Read-Eval-Print Loop). It looks like this:
   ```
   julia>
   ```

5. Enter the package manager mode by pressing `]`. The prompt will change to:
   ```
   (@v1.x) pkg>
   ```
   The `@v1.x` shows your Julia version number.

6. In the package manager, run the `instantiate` command:
   ```
   (@v1.x) pkg> instantiate
   ```
   This command reads the `Project.toml` file in your project directory, creates the `Manifest.toml` file, and installs all the required packages and their dependencies. It might take a few minutes, especially on the first run.

7. Once the instantiation is complete, exit the package manager by pressing the Backspace key. This will return you to the Julia REPL.

8. You can now exit Julia by pressing `ctrl+d` or by typing:
   ```
   julia> exit()
   ```

After completing these steps, your project environment is set up and ready to use. You won't need to repeat this process unless you move the project to a new machine or update the project's dependencies.


## Input Files

- Input files contain the list of edges of the network, one per line.
- Nodes can be specified as integers or strings.
- The program automatically detects the number of columns:
  - If more than two columns are found, the network is treated as weighted.
- File extension (e.g., `.net`) is arbitrary; files are essentially CSV or TSV format.

## Running the Program

Execute the program with the following command:

```
julia --project FitnessCentrality.jl filename
```

Where `filename` is the path to the network file you want to analyze. For example:

```
julia --project FitnessCentrality.jl networks/star-wheel.net
```

### Directed Networks

The program will ask if the network is directed. Type `y` if this is the case.

## Results

The analysis results are placed in the `results/` folder.