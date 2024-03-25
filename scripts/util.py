import os
import pandas as pd
import pickle


def create_folder(path):
        if not os.path.exists(path):
            os.makedirs(path)


# TODO apply this selection at the beginning of each representation
def filter_dict_by_string(dictionary, substring):
    return {key: value for key, value in dictionary.items() if substring in key}


def filter_dict_by_substrings(dictionary, substrings):
    return {key: value for key, value in dictionary.items() if any(substring in key for substring in substrings)}


def set_working_directory(new_directory):
    try:
        os.chdir(new_directory)
        print(f"Working directory changed to: {os.getcwd()}")
    except OSError as e:
        print(f"Error changing working directory: {e}")


def file_exists(file_path):
    try:
        # Attempt to get file information
        file_info = os.stat(file_path)
        return True
    except FileNotFoundError:
        # Handle the case when the file is not found
        return False
    except Exception as e:
        # Handle other exceptions, if any
        print(f"An error occurred: {e}")
        return False


def run_script(script, stdin=None):
    """Returns (stdout, stderr), raises error on non-zero return code"""
    import subprocess
    # Note: by using a list here (['bash', ...]) you avoid quoting issues, as the 
    # arguments are passed in exactly this order (spaces, quotes, and newlines won't
    # cause problems):
    proc = subprocess.Popen(['bash', '-c', script],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        stdin=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    if proc.returncode:
        raise ScriptException(proc.returncode, stdout, stderr, script)
    return stdout, stderr


def run_script_sp(script, stdin=None):
    """Returns (stdout, stderr), raises error on non-zero return code"""
    # Execute the script using os.system
    returncode = os.system(script)
    stdout = None  # Since os.system doesn't capture stdout
    stderr = None  # Since os.system doesn't capture stderr
    if returncode != 0:
        raise ScriptException(returncode, stdout, stderr, script)
    return stdout, stderr


class ScriptException(Exception):
    def __init__(self, returncode, stdout, stderr, script):
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr
        Exception().__init__('Error in script')


def update_pickle_results(simulation_name, output_folder, replica_number, df_toAdd, simulation_prefix, data_prefix):
    
    pickle_file_path = f'{output_folder}/results.pkl'


    # if overResidues is True:
    #     pickle_file_path = f'{output_folder}/{replica_number}/results_overResidue_{replica_number}.pkl'
    # if overTime is True:
    #     pickle_file_path = f'{output_folder}/{replica_number}/results_overTime_{replica_number}.pkl'
    # if overLigand is True:
    #     pickle_file_path = f'{output_folder}/{replica_number}/results_overLigand_{replica_number}.pkl'
    
    try:
        with open(pickle_file_path, 'rb') as file:
            dataframes_dict = pickle.load(file)
    except FileNotFoundError:
        print(f'No output file for {output_folder}. Creating a new data file.')
        dataframes_dict = {}

    dataframes_dict[f'{simulation_name}|{replica_number}|{simulation_prefix}|{data_prefix}'] = df_toAdd

    with open(pickle_file_path, 'wb') as file:
        pickle.dump(dataframes_dict, file)


# def update_results_df(output_folder, replica_number, df_toAdd, column_prefix, overResidues = True, overTime = False, overLigand = False):
#     df_toAdd = df_toAdd.add_prefix(f'{column_prefix}:', axis='columns')
#     if overResidues is True:
#         csv_file_path = f'{output_folder}/{replica_number}/results_overResidue_{replica_number}.csv'
#     if overTime is True:
#         csv_file_path = f'{output_folder}/{replica_number}/results_overTime_{replica_number}.csv'
#     if overLigand is True:
#         csv_file_path = f'{output_folder}/{replica_number}/results_overLigand_{replica_number}.csv'
#     try:
#         # Read the existing CSV file into a DataFrame
#         df = pd.read_csv(csv_file_path, index_col=0)
#     except FileNotFoundError:
#         # Create an empty DataFrame with new columns and data
#         df = pd.DataFrame()
    
#     dfs_to_concat = [df]  # Initialize with the existing DataFrame
#     for column, data in df_toAdd.items():
#         if column in df.columns:
#             # Column exists, update its values
#             df[column] = data
#         else:
#             # Column doesn't exist, add it to the list of DataFrames
#             dfs_to_concat.append(pd.DataFrame({column: data}))

#     # Concatenate DataFrames along columns
#     df = pd.concat(dfs_to_concat, axis=1)

#     if overTime is True:
#         df.sort_index(inplace=True)

#     # Save the DataFrame to the CSV file
#     df.to_csv(csv_file_path, index=True)


def st_get_simulation_type(column_list):
    simulation_name, simulation_type, data_type = [], [], []
    for column in column_list:
        name, simulation, data = column.split(':')
        simulation_name.append(name)
        simulation_type.append(simulation)
        data_type.append(data)
    return set(simulation_name), set(simulation_type), set(data_type)


def st_get_simulation_info(selected_simulations_dict):
    simulation_names, simulation_replicas, simulation_subsets, data_type = set(), set(), set(), set()
    
    for k in selected_simulations_dict.keys():
        name, replica, subset, data = k.split('|')
        
        simulation_names.add(name)
        simulation_replicas.add(replica)
        simulation_subsets.add(subset)
        data_type.add(data)

    return sorted(list(simulation_names)), sorted(list(simulation_replicas)), sorted(list(simulation_subsets)), sorted(list(data_type))
