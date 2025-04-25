import sys
import os
import glob
import csv
import matplotlib.pyplot as plt
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QListWidget, QPushButton, QCheckBox, QLabel, QInputDialog,
    QScrollArea, QSlider, QMessageBox
)
from PyQt5.QtCore import Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

# Enable dark mode for matplotlib plots
plt.style.use("seaborn-v0_8-dark-palette")

class PlotCanvas(FigureCanvas):
    def __init__(self, parent=None):
        fig, self.ax = plt.subplots()
        super().__init__(fig)
        self.setParent(parent)
        self.title_text = None
        self.x_label_text = None
        self.y_label_text = None

        # Connect pick events to handle clicks on text elements
        self.figure.canvas.mpl_connect('pick_event', self.on_pick)

    def plot(self, x_axis, y_axes, x_title, y_titles, title):
        self.ax.clear()
        for y_axis, y_title in zip(y_axes, y_titles):
            self.ax.plot(x_axis, y_axis, label=y_title)
        # Set labels and make them pickable
        self.x_label_text = self.ax.set_xlabel(x_title, picker=True)
        self.y_label_text = self.ax.set_ylabel("Variables", picker=True)
        self.title_text = self.ax.set_title(title, picker=True)
        self.ax.legend()
        self.draw()

    def on_pick(self, event):
        artist = event.artist
        if artist == self.title_text:
            self.edit_title()
        elif artist == self.x_label_text:
            self.edit_x_label()
        elif artist == self.y_label_text:
            self.edit_y_label()

    def edit_title(self):
        text, ok = QInputDialog.getText(self.parent(), "Edit Title", "Enter new title:")
        if ok and text:
            self.title_text.set_text(text)
            self.draw()

    def edit_x_label(self):
        text, ok = QInputDialog.getText(self.parent(), "Edit X-Axis Label", "Enter new x-axis label:")
        if ok and text:
            self.x_label_text.set_text(text)
            self.draw()

    def edit_y_label(self):
        text, ok = QInputDialog.getText(self.parent(), "Edit Y-Axis Label", "Enter new y-axis label:")
        if ok and text:
            self.y_label_text.set_text(text)
            self.draw()

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("CSV Plotter - Dark Mode")
        self.setGeometry(100, 100, 1000, 600)
        self.setStyleSheet("background-color: #2e2e2e; color: white;")

        # Dictionary to store full paths with filenames as keys
        self.file_paths = {}
        self.current_file_name = ""
        self.current_spectrum = []
        self.current_time_list = []
        self.selected_steps = 1  # Default value for the slider
        self.current_vars = []   # Store variables selected for saving

        # Main widget and layout
        main_widget = QWidget(self)
        self.setCentralWidget(main_widget)
        main_layout = QHBoxLayout()
        main_widget.setLayout(main_layout)

        # Left panel for file list and selection feedback
        left_layout = QVBoxLayout()
        self.file_list = QListWidget()
        self.file_list.setStyleSheet(
            "QListWidget { background-color: #2e2e2e; color: white; }"
        )
        self.load_csv_files()
        self.file_list.itemClicked.connect(self.load_and_display_vars)
        left_layout.addWidget(self.file_list)

        # Variable selection message label
        self.message_label = QLabel("")
        self.message_label.setStyleSheet("color: red;")
        left_layout.addWidget(self.message_label)

        # Number of variables label
        self.num_vars_label = QLabel("Number of variables (excluding Time): 0")
        self.num_vars_label.setStyleSheet("color: white;")
        left_layout.addWidget(self.num_vars_label)

        # Checkbox container for variables
        self.vars_container = QScrollArea()
        self.vars_container.setWidgetResizable(True)
        self.vars_widget = QWidget()
        self.vars_layout = QVBoxLayout()
        self.vars_widget.setLayout(self.vars_layout)
        self.vars_container.setWidget(self.vars_widget)
        left_layout.addWidget(self.vars_container)

        # Unselect All button
        self.unselect_all_button = QPushButton("Unselect All")
        self.unselect_all_button.clicked.connect(self.unselect_all_vars)
        self.unselect_all_button.setStyleSheet(
            "QPushButton { background-color: #3a3a3a; color: white; border: none; }"
            "QPushButton:hover { background-color: #505050; }"
        )
        left_layout.addWidget(self.unselect_all_button)

        # Slider for time steps
        self.slider_label = QLabel("Number of Time Steps: 1")
        self.slider_label.setStyleSheet("color: white;")
        self.time_slider = QSlider(Qt.Horizontal)
        self.time_slider.setMinimum(1)
        self.time_slider.valueChanged.connect(self.update_slider_label)
        self.time_slider.valueChanged.connect(self.plot_selected_vars)  # Auto-update plot on slider change
        left_layout.addWidget(self.slider_label)
        left_layout.addWidget(self.time_slider)

        # Checkbox for ignoring zero values
        self.ignore_zero_checkbox = QCheckBox("Ignore zero values")
        self.ignore_zero_checkbox.setStyleSheet("QCheckBox { color: white; }")
        self.ignore_zero_checkbox.stateChanged.connect(self.load_and_display_vars)
        left_layout.addWidget(self.ignore_zero_checkbox)

        # Checkbox for ignoring duplicate variables
        self.ignore_duplicates_checkbox = QCheckBox("Ignore duplicates")
        self.ignore_duplicates_checkbox.setStyleSheet("QCheckBox { color: white; }")
        self.ignore_duplicates_checkbox.stateChanged.connect(self.load_and_display_vars)
        left_layout.addWidget(self.ignore_duplicates_checkbox)

        # Checkbox for default naming
        self.default_name_checkbox = QCheckBox("Use default name")
        self.default_name_checkbox.setStyleSheet("QCheckBox { color: white; }")
        left_layout.addWidget(self.default_name_checkbox)

        # Save Plot button
        self.save_button = QPushButton("Save Plot")
        self.save_button.clicked.connect(self.save_plot)
        self.save_button.setStyleSheet(
            "QPushButton { background-color: #3a3a3a; color: white; border: none; }"
            "QPushButton:hover { background-color: #505050; }"
        )
        left_layout.addWidget(self.save_button)

        # Canvas for displaying plots
        self.canvas = PlotCanvas(self)
        right_layout = QVBoxLayout()
        right_layout.addWidget(self.canvas)

        # Adding layouts to main layout
        main_layout.addLayout(left_layout, 1)
        main_layout.addLayout(right_layout, 3)

    def load_csv_files(self):
        # Load CSV files sorted by modification time
        csv_files = glob.glob(os.path.join(os.getcwd(), "data", "csv", "*.csv"))
        csv_files.sort(key=os.path.getmtime, reverse=True)
        for file_path in csv_files:
            file_name = os.path.splitext(os.path.basename(file_path))[0]  # Extract file name without extension
            self.file_list.addItem(file_name)  # Add file name to list
            self.file_paths[file_name] = file_path  # Store full path in dictionary

    def read_csv(self, filename):
        file_path = self.file_paths[filename]  # Retrieve full path from dictionary
        with open(file_path, "r") as f:
            reader = csv.reader(f)
            spectrum = [list(map(float, row)) for row in reader]
            time_list = spectrum.pop(-1)
            return self.transpose(spectrum), time_list

    def transpose(self, matrix):
        return [[matrix[j][i] for j in range(len(matrix))] for i in range(len(matrix[0]))]

    def load_and_display_vars(self):
        # Check if a file is selected
        selected_item = self.file_list.currentItem()
        if selected_item:
            self.current_file_name = selected_item.text()
            self.current_spectrum, self.current_time_list = self.read_csv(self.current_file_name)

            # Set the slider range to the number of time steps
            num_timesteps = len(self.current_time_list)
            self.time_slider.setMaximum(num_timesteps)
            self.time_slider.setValue(num_timesteps)  # Default to max time steps
        elif not self.current_file_name:
            return  # No file selected, do nothing

        # Clear any existing checkboxes in the vars_layout
        for i in reversed(range(self.vars_layout.count())):
            widget = self.vars_layout.itemAt(i).widget()
            if widget:
                widget.setParent(None)

        # Reset current_vars
        self.current_vars = []

        # Add "Time" checkbox
        time_checkbox = QCheckBox("Time")
        time_checkbox.setStyleSheet("QCheckBox { color: white; }")
        time_checkbox.setChecked(False)  # Ensure it's unchecked
        time_checkbox.stateChanged.connect(self.plot_selected_vars)  # Auto-update plot on checkbox change
        self.vars_layout.addWidget(time_checkbox)

        num_vars_displayed = 0  # Counter for number of variables displayed (excluding Time)
        seen_vars = set()  # Set to store seen variables if ignoring duplicates

        # Add checkboxes for each variable
        for i, var_data in enumerate(self.current_spectrum):
            # If 'Ignore zero values' is checked, skip variables that are all zero
            if self.ignore_zero_checkbox.isChecked():
                if all(v == 0 for v in var_data):
                    continue

            # If 'Ignore duplicates' is checked, skip duplicate variables
            if self.ignore_duplicates_checkbox.isChecked():
                # Quantize data to 4 decimal places to avoid floating point precision issues
                var_data_quantized = tuple(round(v, 4) for v in var_data)
                if var_data_quantized in seen_vars:
                    continue
                else:
                    seen_vars.add(var_data_quantized)

            checkbox = QCheckBox(f"x{i + 1}")
            checkbox.setStyleSheet("QCheckBox { color: white; }")
            checkbox.setChecked(False)  # Ensure it's unchecked
            checkbox.stateChanged.connect(self.plot_selected_vars)  # Auto-update plot on checkbox change
            self.vars_layout.addWidget(checkbox)
            num_vars_displayed += 1

        # Update the number of variables label
        self.num_vars_label.setText(f"Number of variables (excluding Time): {num_vars_displayed}")

    def unselect_all_vars(self):
        for cb in self.vars_widget.findChildren(QCheckBox):
            cb.setChecked(False)

    def update_slider_label(self, value):
        # Update the label and set the selected time steps
        self.slider_label.setText(f"Number of Time Steps: {value}")
        self.selected_steps = value

    def plot_selected_vars(self):
        # Get selected variables
        selected_vars = []
        for cb in self.vars_widget.findChildren(QCheckBox):
            if cb.isChecked():
                if cb.text() == "Time":
                    selected_vars.append("Time")
                else:
                    selected_vars.append(int(cb.text().replace("x", "")) - 1)

        # Validate the selection before plotting
        if "Time" in selected_vars:
            if len(selected_vars) < 2:
                self.message_label.setText("Select at least one variable in addition to Time.")
                return
            self.message_label.setText("")  # Clear message if selection is valid
            x_axis = self.current_time_list[:self.selected_steps]
            y_axes = [self.current_spectrum[var][:self.selected_steps] for var in selected_vars if var != "Time"]
            x_title = "Time"
            y_titles = [f"x{var + 1}" for var in selected_vars if var != "Time"]
        else:
            if len(selected_vars) != 2:
                self.message_label.setText("Select exactly two variables for plotting.")
                return
            self.message_label.setText("")  # Clear message if selection is valid
            var1, var2 = selected_vars
            x_axis = self.current_spectrum[var1][:self.selected_steps]
            y_axes = [self.current_spectrum[var2][:self.selected_steps]]
            x_title = f"x{var1 + 1}"
            y_titles = [f"x{var2 + 1}"]

        title = f"Equations of Motion for {self.current_file_name}"

        # Plotting
        self.canvas.plot(x_axis, y_axes, x_title, y_titles, title)

        # Store variables for saving with names that match the plot titles
        self.current_vars = [x_title] + y_titles

    def save_plot(self):
        if not self.current_vars:
            QMessageBox.warning(self, "Warning", "Please select variables to plot before saving.")
            return

        var1, *vars_rest = self.current_vars
        if self.default_name_checkbox.isChecked():
            save_path = os.path.join(
                os.getcwd(), "data", "plots",
                f"{self.current_file_name}_{var1}_{'_'.join(vars_rest)}.png"
            )
        else:
            custom_name, ok = QInputDialog.getText(self, "Save As", "Enter custom name for the plot:")
            if not ok or not custom_name:
                return
            save_path = os.path.join(os.getcwd(), "data", "plots", f"{custom_name}.png")

        # Use a dark background for the saved plot
        self.canvas.figure.savefig(save_path, facecolor="#2e2e2e")
        QMessageBox.information(self, "Plot Saved", f"Plot saved as: {save_path}")

# Run the application
if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
