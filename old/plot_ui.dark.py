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
plt.style.use("seaborn-v0_8-white")

# Set default text and axis colors for dark background
plt.rcParams['text.color'] = 'white'
plt.rcParams['axes.labelcolor'] = 'white'
plt.rcParams['xtick.color'] = 'white'
plt.rcParams['ytick.color'] = 'white'
plt.rcParams['axes.edgecolor'] = 'white'

class PlotCanvas(FigureCanvas):
    def __init__(self, parent=None):
        fig, self.ax = plt.subplots()
        # Set figure and axes facecolor to match dark mode
        fig.patch.set_facecolor('#2e2e2e')
        self.ax.set_facecolor('#2e2e2e')
        super().__init__(fig)
        self.setParent(parent)
        self.title_text = None
        self.x_label_text = None
        self.y_label_text = None

        # Connect pick events to handle clicks on text elements
        self.figure.canvas.mpl_connect('pick_event', self.on_pick)

    def plot(self, x_axis, y_axes, x_title, y_titles, title):
        self.ax.clear()
        # Ensure facecolors remain dark after clear
        self.ax.set_facecolor('#2e2e2e')
        for y_axis, y_title in zip(y_axes, y_titles):
            self.ax.plot(x_axis, y_axis, label=y_title)
        # Set labels with explicit color and make them pickable
        self.x_label_text = self.ax.set_xlabel(x_title, picker=True, color='white')
        self.y_label_text = self.ax.set_ylabel("Variables", picker=True, color='white')
        self.title_text = self.ax.set_title(title, picker=True, color='white')
        # Update tick and spine colors
        self.ax.tick_params(axis='x', colors='white')
        self.ax.tick_params(axis='y', colors='white')
        for spine in self.ax.spines.values():
            spine.set_color('white')
        # Adjust legend appearance
        legend = self.ax.legend(facecolor='#3a3a3a', edgecolor='white')
        for text in legend.get_texts():
            text.set_color('white')
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
            self.title_text.set_color('white')
            self.draw()

    def edit_x_label(self):
        text, ok = QInputDialog.getText(self.parent(), "Edit X-Axis Label", "Enter new x-axis label:")
        if ok and text:
            self.x_label_text.set_text(text)
            self.x_label_text.set_color('white')
            self.draw()

    def edit_y_label(self):
        text, ok = QInputDialog.getText(self.parent(), "Edit Y-Axis Label", "Enter new y-axis label:")
        if ok and text:
            self.y_label_text.set_text(text)
            self.y_label_text.set_color('white')
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
        self.traced_flag = False # Will be set when a "Traced" file is loaded

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
        self.time_slider.valueChanged.connect(self.plot_selected_vars)
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
            file_name = os.path.splitext(os.path.basename(file_path))[0]
            self.file_list.addItem(file_name)
            self.file_paths[file_name] = file_path

    def read_csv(self, filename):
        file_path = self.file_paths[filename]
        with open(file_path, "r") as f:
            reader = csv.reader(f)
            spectrum = [list(map(float, row)) for row in reader]
            time_list = spectrum.pop(-1)
            return self.transpose(spectrum), time_list

    def transpose(self, matrix):
        return [[matrix[j][i] for j in range(len(matrix))] for i in range(len(matrix[0]))]

    def load_and_display_vars(self):
        selected_item = self.file_list.currentItem()
        if selected_item:
            self.current_file_name = selected_item.text()
            # detect traced files
            self.traced_flag = "traced" in self.current_file_name.lower()

            self.current_spectrum, self.current_time_list = self.read_csv(self.current_file_name)
            num_timesteps = len(self.current_time_list)
            self.time_slider.setMaximum(num_timesteps)
            self.time_slider.setValue(num_timesteps)
        elif not self.current_file_name:
            return

        # clear existing checkboxes
        for i in reversed(range(self.vars_layout.count())):
            widget = self.vars_layout.itemAt(i).widget()
            if widget:
                widget.setParent(None)

        self.current_vars = []
        # Time checkbox
        time_checkbox = QCheckBox("Time")
        time_checkbox.setStyleSheet("QCheckBox { color: white; }")
        time_checkbox.stateChanged.connect(self.plot_selected_vars)
        self.vars_layout.addWidget(time_checkbox)

        num_vars_displayed = 0
        seen_vars = set()

        for i, var_data in enumerate(self.current_spectrum):
            if self.ignore_zero_checkbox.isChecked() and all(v == 0 for v in var_data):
                continue
            if self.ignore_duplicates_checkbox.isChecked():
                quant = tuple(round(v, 4) for v in var_data)
                if quant in seen_vars:
                    continue
                seen_vars.add(quant)

            # custom labeling for traced files
            if self.traced_flag and i < 3:
                labels = ["<<Tr(XX)>>/N", "<<Tr(PP)>>/N", "<<Tr(XP)>>/N"]
                cb_label = labels[i]
            else:
                cb_label = f"x{i+1}"

            checkbox = QCheckBox(cb_label)
            checkbox.setStyleSheet("QCheckBox { color: white; }")
            checkbox.stateChanged.connect(self.plot_selected_vars)
            self.vars_layout.addWidget(checkbox)
            num_vars_displayed += 1

        self.num_vars_label.setText(f"Number of variables (excluding Time): {num_vars_displayed}")

    def unselect_all_vars(self):
        for cb in self.vars_widget.findChildren(QCheckBox):
            cb.setChecked(False)

    def update_slider_label(self, value):
        self.slider_label.setText(f"Number of Time Steps: {value}")
        self.selected_steps = value

    def plot_selected_vars(self):
        selected = []
        for cb in self.vars_widget.findChildren(QCheckBox):
            if cb.isChecked():
                if cb.text() == "Time":
                    selected.append("Time")
                else:
                    try:
                        selected.append(int(cb.text().replace("x", "")) - 1)
                    except ValueError:
                        selected.append(cb.text())

        if "Time" in selected:
            if len(selected) < 2:
                self.message_label.setText("Select at least one variable in addition to Time.")
                return
            self.message_label.setText("")
            x_axis = self.current_time_list[:self.selected_steps]
            y_axes, y_titles = [], []
            for var in selected:
                if var == "Time":
                    continue
                if isinstance(var, int):
                    y_axes.append(self.current_spectrum[var][:self.selected_steps])
                    y_titles.append(f"x{var+1}")
                else:
                    idx_map = {"<<Tr(XX)>>/N":0, "<<Tr(PP)>>/N":1, "<<Tr(XP)>>/N":2}
                    idx = idx_map[var]
                    y_axes.append(self.current_spectrum[idx][:self.selected_steps])
                    y_titles.append(var)
            x_title = "Time"

        else:
            if len(selected) != 2:
                self.message_label.setText("Select exactly two variables for plotting.")
                return
            self.message_label.setText("")
            var1, var2 = selected
            x_axis = self.current_spectrum[var1][:self.selected_steps]
            y_axes = [self.current_spectrum[var2][:self.selected_steps]]
            x_title = f"x{var1+1}"
            y_titles = [f"x{var2+1}"]

        # dimension mismatch check
        for y in y_axes:
            if len(y) != len(x_axis):
                self.message_label.setText("Error: x and y data have different lengths.")
                return

        title = f"Equations of Motion for {self.current_file_name}"
        self.canvas.plot(x_axis, y_axes, x_title, y_titles, title)
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

        # Save with dark facecolor but ensure text remains visible
        self.canvas.figure.savefig(save_path, facecolor="#2e2e2e", edgecolor="#2e2e2e")
        QMessageBox.information(self, "Plot Saved", f"Plot saved as: {save_path}")

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
