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
        fig.patch.set_facecolor('#2e2e2e')
        self.ax.set_facecolor('#2e2e2e')
        super().__init__(fig)
        self.setParent(parent)
        self.title_text = None
        self.x_label_text = None
        self.y_label_text = None
        self.figure.canvas.mpl_connect('pick_event', self.on_pick)

    def plot(self, x_axis, y_axes, x_title, y_titles, title):
        self.ax.clear()
        self.ax.set_facecolor('#2e2e2e')
        for y_axis, y_title in zip(y_axes, y_titles):
            self.ax.plot(x_axis, y_axis, label=y_title)
        self.x_label_text = self.ax.set_xlabel(x_title, picker=True, color='white')
        self.y_label_text = self.ax.set_ylabel("Variables", picker=True, color='white')
        self.title_text = self.ax.set_title(title, picker=True, color='white')
        self.ax.tick_params(axis='x', colors='white')
        self.ax.tick_params(axis='y', colors='white')
        for spine in self.ax.spines.values():
            spine.set_color('white')
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
        self.file_paths = {}
        self.current_file_name = ""
        self.current_spectrum = []
        self.current_time_list = []
        self.selected_steps = 1
        self.traced_flag = False

        main_widget = QWidget(self)
        self.setCentralWidget(main_widget)
        main_layout = QHBoxLayout()
        main_widget.setLayout(main_layout)

        # Left panel
        left_layout = QVBoxLayout()
        self.file_list = QListWidget()
        self.file_list.setStyleSheet("QListWidget { background-color: #2e2e2e; color: white; }")
        self.load_csv_files()
        self.file_list.itemClicked.connect(self.load_and_display_vars)
        left_layout.addWidget(self.file_list)

        self.message_label = QLabel("")
        self.message_label.setStyleSheet("color: red;")
        left_layout.addWidget(self.message_label)

        self.num_vars_label = QLabel("Number of variables (excluding Time): 0")
        self.num_vars_label.setStyleSheet("color: white;")
        left_layout.addWidget(self.num_vars_label)

        self.vars_container = QScrollArea()
        self.vars_container.setWidgetResizable(True)
        self.vars_widget = QWidget()
        self.vars_layout = QVBoxLayout()
        self.vars_widget.setLayout(self.vars_layout)
        self.vars_container.setWidget(self.vars_widget)
        left_layout.addWidget(self.vars_container)

        btn_layout = QHBoxLayout()
        self.ignore_zero_checkbox = QCheckBox("Ignore zero values")
        self.ignore_zero_checkbox.setStyleSheet("QCheckBox { color: white; }")
        self.ignore_zero_checkbox.stateChanged.connect(self.load_and_display_vars)
        btn_layout.addWidget(self.ignore_zero_checkbox)

        self.ignore_duplicates_checkbox = QCheckBox("Ignore duplicates")
        self.ignore_duplicates_checkbox.setStyleSheet("QCheckBox { color: white; }")
        self.ignore_duplicates_checkbox.stateChanged.connect(self.load_and_display_vars)
        btn_layout.addWidget(self.ignore_duplicates_checkbox)
        left_layout.addLayout(btn_layout)

        self.time_slider = QSlider(Qt.Horizontal)
        self.time_slider.setMinimum(1)
        self.time_slider.valueChanged.connect(self.update_slider_label)
        self.time_slider.valueChanged.connect(self.plot_selected_vars)
        left_layout.addWidget(self.time_slider)

        self.slider_label = QLabel("Number of Time Steps: 1")
        self.slider_label.setStyleSheet("color: white;")
        left_layout.addWidget(self.slider_label)

        self.unselect_all_button = QPushButton("Unselect All")
        self.unselect_all_button.clicked.connect(self.unselect_all_vars)
        self.unselect_all_button.setStyleSheet(
            "QPushButton { background-color: #3a3a3a; color: white; border: none; }"
            "QPushButton:hover { background-color: #505050; }"
        )
        left_layout.addWidget(self.unselect_all_button)

        self.default_name_checkbox = QCheckBox("Use default name")
        self.default_name_checkbox.setStyleSheet("QCheckBox { color: white; }")
        left_layout.addWidget(self.default_name_checkbox)

        self.save_button = QPushButton("Save Plot")
        self.save_button.clicked.connect(self.save_plot)
        self.save_button.setStyleSheet(
            "QPushButton { background-color: #3a3a3a; color: white; border: none; }"
            "QPushButton:hover { background-color: #505050; }"
        )
        left_layout.addWidget(self.save_button)

        # Right panel
        self.canvas = PlotCanvas(self)
        right_layout = QVBoxLayout()
        right_layout.addWidget(self.canvas)

        main_layout.addLayout(left_layout, 1)
        main_layout.addLayout(right_layout, 3)

    def load_csv_files(self):
        csv_dir = os.path.join(os.getcwd(), "data", "csv")
        files = glob.glob(os.path.join(csv_dir, "*.csv"))
        files.sort(key=os.path.getmtime, reverse=True)
        for fp in files:
            name = os.path.splitext(os.path.basename(fp))[0]
            self.file_list.addItem(name)
            self.file_paths[name] = fp

    def read_csv(self, name):
        with open(self.file_paths[name], 'r') as f:
            reader = csv.reader(f)
            data = [list(map(float, row)) for row in reader]
        time = data.pop(-1)
        return self.transpose(data), time

    def transpose(self, M):
        return list(map(list, zip(*M)))

    def load_and_display_vars(self):
        item = self.file_list.currentItem()
        if not item: return
        self.current_file_name = item.text()
        self.traced_flag = "traced" in self.current_file_name.lower()
        self.current_spectrum, self.current_time_list = self.read_csv(self.current_file_name)
        self.time_slider.setMaximum(len(self.current_time_list))
        self.time_slider.setValue(len(self.current_time_list))

        # rebuild variable checkboxes
        for i in reversed(range(self.vars_layout.count())):
            w = self.vars_layout.itemAt(i).widget()
            if w: w.setParent(None)

        # "Time" checkbox
        cb = QCheckBox("Time")
        cb.setStyleSheet("QCheckBox { color: white; }")
        cb.stateChanged.connect(self.plot_selected_vars)
        self.vars_layout.addWidget(cb)

        seen = set()
        count = 0
        for idx, series in enumerate(self.current_spectrum):
            if self.ignore_zero_checkbox.isChecked() and all(v==0 for v in series):
                continue
            if self.ignore_duplicates_checkbox.isChecked():
                key = tuple(round(v,4) for v in series)
                if key in seen:
                    continue
                seen.add(key)

            label = f"x{idx+1}"
            if self.traced_flag and idx<3:
                label = ["<<Tr(XX)>>/N","<<Tr(PP)>>/N","<<Tr(XP)>>/N"][idx]

            cb = QCheckBox(label)
            cb.setStyleSheet("QCheckBox { color: white; }")
            cb.stateChanged.connect(self.plot_selected_vars)
            self.vars_layout.addWidget(cb)
            count += 1

        self.num_vars_label.setText(f"Number of variables (excluding Time): {count}")

    def unselect_all_vars(self):
        for cb in self.vars_widget.findChildren(QCheckBox):
            cb.setChecked(False)

    def update_slider_label(self, v):
        self.selected_steps = v
        self.slider_label.setText(f"Number of Time Steps: {v}")

    def plot_selected_vars(self):
        # gather selections
        sels = []
        for cb in self.vars_widget.findChildren(QCheckBox):
            if cb.isChecked():
                sels.append(cb.text())
        if "Time" in sels:
            if len(sels)<2:
                self.message_label.setText("Select at least one variable in addition to Time.")
                return
            self.message_label.setText("")
            xdata = self.current_time_list[:self.selected_steps]
            ydatas, ylabels = [], []
            for text in sels:
                if text=="Time": continue
                if text.startswith("x"):
                    i = int(text[1:])-1
                else:
                    i = {"<<Tr(XX)>>/N":0,"<<Tr(PP)>>/N":1,"<<Tr(XP)>>/N":2}[text]
                ydatas.append(self.current_spectrum[i][:self.selected_steps])
                ylabels.append(text)
            xlabel = "Time"
        else:
            if len(sels)!=2:
                self.message_label.setText("Select exactly two variables for plotting.")
                return
            self.message_label.setText("")
            i1 = int(sels[0].lstrip("x"))-1
            i2 = int(sels[1].lstrip("x"))-1
            xdata = self.current_spectrum[i1][:self.selected_steps]
            ydatas = [self.current_spectrum[i2][:self.selected_steps]]
            xlabel = sels[0]
            ylabels = [sels[1]]

        title = f"Equations of Motion for {self.current_file_name}"
        self.canvas.plot(xdata, ydatas, xlabel, ylabels, title)

    def save_plot(self):
        # must have plotted at least once
        if not hasattr(self.canvas, 'title_text') or self.canvas.title_text is None:
            QMessageBox.warning(self, "Warning", "Plot something first.")
            return

        # build filename
        if self.default_name_checkbox.isChecked():
            fname = f"{self.current_file_name}_{self.canvas.x_label_text.get_text()}_" + \
                    "_".join(t.get_text() for t in self.canvas.ax.get_legend().get_texts()) + ".png"
        else:
            name, ok = QInputDialog.getText(self, "Save As", "Enter filename:")
            if not ok or not name:
                return
            fname = f"{name}.png"

        outdir = os.path.join(os.getcwd(), "data", "plots")
        os.makedirs(outdir, exist_ok=True)
        outpath = os.path.join(outdir, fname)

        # pull the current data out of the canvas
        xdata = self.canvas.ax.lines[0].get_xdata()
        ydatas = [line.get_ydata() for line in self.canvas.ax.lines]
        # read labels from canvas
        title = self.canvas.title_text.get_text()
        xlabel = self.canvas.x_label_text.get_text()
        ylabel = self.canvas.y_label_text.get_text()
        ylabels = [t.get_text() for t in self.canvas.ax.get_legend().get_texts()]

        # redraw in light style and save
        with plt.style.context("seaborn-v0_8-white"):
            fig, ax = plt.subplots()
            for ydata, label in zip(ydatas, ylabels):
                ax.plot(xdata, ydata, label=label)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.set_title(title)
            ax.legend()
            fig.savefig(outpath, dpi=150)
            plt.close(fig)

        QMessageBox.information(self, "Saved", f"Plot saved to:\n{outpath}")

if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = MainWindow()
    w.show()
    sys.exit(app.exec_())
