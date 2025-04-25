import matplotlib.pyplot as plt
import csv
import os
import glob
import sys
'''
This script is used to plot the latest csv file in the data/csv directory. Naming scheme is used for undestanding what to plot.
- If the csv file is EOM, then it plots the equations of motion, it prompts user for which variables to plot.
  - It asks user for x1, x2, x3, ...., t. t is the latest line of the csv file, and the other variables are the columns of the csv file, excluded the last line.
- If the csv file is LYAP, then it plots the lyapunov exponents vs the energy.
-  Then it saves the plot in the data/plots directory and shows it.
'''

def transpose(matrix):
  '''
  Returns the transpose of a matrix, a matrix is a list of lists.
  '''
  return [[matrix[j][i] for j in range(len(matrix))] for i in range(len(matrix[0]))]

def get_latest_csv():
  '''
  Return the latest csv file in the data/csv directory.
  '''
  csv_files = glob.glob(os.path.join(os.getcwd(), 'data', 'csv', '*.csv'))
  # Sorting the files by last editing time then getting the latest editted csv file.
  csv_files.sort(key=os.path.getmtime)
  return csv_files[-1]

def read_csv(filename):
  '''
  Reads the csv file and returns the spectrum and time_list.
  '''
  with open(filename, 'r') as f:
    reader = csv.reader(f)
    spectrum = [list(map(float, row)) for row in reader]
    time_list = spectrum.pop(-1)
    return transpose(spectrum), time_list

def plot_vars(x_axis, y_axis, x_title, y_title, title, filename):
  # plt.style.use('seaborn')
  plt.plot(x_axis, y_axis)
  plt.xlabel(x_title)
  plt.ylabel(y_title)
  plt.title(title)
  plt.savefig(os.path.join(os.getcwd(), 'data', 'plots', f'{filename}_{x_title}_{y_title}.png'))
  plt.show()

def plot_eom(spectrum, time_list, csv_file_full_path):
  '''
  Plots the equations of motion.
  '''
  filename = os.path.basename(csv_file_full_path).split('.')[0]
  title = f'Equations of Motion for {filename}'
  x_axis = None
  y_axis = None
  x_title = None
  y_title = None
  vars = input('Which variables do you want to plot? (x0, x1, x2, x3, ..., xn): ')
  try:
    vars = [int (a) for a in vars.replace('x', '').strip().split(',')]
    if vars == []:
      print('No variables selected.')
      plot_eom(spectrum, time_list, csv_file_full_path)
    elif len(vars) > 2:
      print('Cannot plot more than 2 variables.')
      plot_eom(spectrum, time_list, csv_file_full_path)
    elif len(vars) == 1:
      print('Need to select 2 variables.')
      plot_eom(spectrum, time_list, csv_file_full_path)
    elif vars[0] > len(spectrum):
      print('Variable does not exist.')
      plot_eom(spectrum, time_list, csv_file_full_path)
    elif vars[1] > len(spectrum):
      print('Variable does not exist.')
      plot_eom(spectrum, time_list, csv_file_full_path)
    elif 0 in vars:
      x_axis = time_list
      vars.remove(0)
      y_axis = spectrum[vars[0] - 1]
      x_title = 'Time'
      y_title = f'x{vars[0]}'
    else:
      x_axis = spectrum[vars[0] - 1]
      y_axis = spectrum[vars[1] - 1]
      x_title = f'x{vars[0]}'
      y_title = f'x{vars[1]}'
    plot_vars(x_axis, y_axis, x_title, y_title, title, filename)
  except:
    print('Invalid input.')
    plot_eom(spectrum, time_list, csv_file_full_path)

if __name__ == '__main__':
  csv_file = get_latest_csv()
  spectrum, time_list = read_csv(csv_file)
  if 'EOM' in csv_file:
    plot_eom(spectrum, time_list, csv_file)
  # elif 'LYAP' in csv_file:
  #   plot.plot_lyap(csv_file)
  else:
    print('Unknown file type.')
