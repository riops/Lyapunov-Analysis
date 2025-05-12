import sys
import os
import glob
import csv
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QListWidget, QPushButton, QCheckBox, QLabel, QInputDialog,
    QScrollArea, QSlider, QMessageBox
)
from PyQt5.QtCore import Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

# Dark palette + MathText via $...$
plt.style.use("seaborn-v0_8-dark-palette")
mpl.rcParams['text.usetex'] = False

# ─── Index‐inversion utilities ─────────────────────────────────────────────

def reverse_indexX(idx: int, N: int):
    block = N * N
    group, offset = divmod(idx, block)
    a = group // 2 + 1
    i = group % 2  + 1
    l = math.isqrt(offset)
    m = offset - (l*l + l)
    return a, i, l, m

def reverse_indexXX(idx: int, N: int):
    max_index = 4 * (N * N)
    i1, i2 = divmod(idx, max_index)
    a1, i1_, l1, m1 = reverse_indexX(i1, N)
    a2, i2_, l2, m2 = reverse_indexX(i2, N)
    return a1, i1_, l1, m1, a2, i2_, l2, m2

def generate_var_label(idx: int, N: int):
    num_per_block = 16 * (N**4)
    block_idx      = idx // num_per_block       # 0=XX,1=PP,2=XP
    inner_idx      = idx % num_per_block
    names          = ['XX','PP','XP']
    sym_map        = {'XX':('X','X'), 'PP':('P','P'), 'XP':('X','P')}
    blk = names[block_idx]
    S1, S2 = sym_map[blk]
    a1,i1,l1,m1,a2,i2,l2,m2 = reverse_indexXX(inner_idx, N)
    return rf"$\langle {S1}_{{{i1},{a1}}}^{{{l1},{m1}}}\,{S2}_{{{i2},{a2}}}^{{{l2},{m2}}}\rangle$"

# ─── Plot canvas ────────────────────────────────────────────────────────────

class PlotCanvas(FigureCanvas):
    def __init__(self, parent=None):
        fig, self.ax = plt.subplots()
        super().__init__(fig)
        self.setParent(parent)
        self.title_text = None
        self.x_label_text = None
        self.y_label_text = None
        self.figure.canvas.mpl_connect('pick_event', self.on_pick)

    def plot(self, x_axis, y_axes, x_title, y_titles, title):
        self.ax.clear()
        for y_axis, y_label in zip(y_axes, y_titles):
            self.ax.plot(x_axis, y_axis, label=y_label)
        self.x_label_text = self.ax.set_xlabel(x_title, picker=True)
        self.y_label_text = self.ax.set_ylabel("Variables", picker=True)
        self.title_text   = self.ax.set_title(title, picker=True)
        self.ax.legend()
        self.draw()

    def on_pick(self, event):
        art = event.artist
        if   art == self.title_text:   self.edit_title()
        elif art == self.x_label_text: self.edit_x_label()
        elif art == self.y_label_text: self.edit_y_label()

    def edit_title(self):
        t, ok = QInputDialog.getText(self.parent(), "Edit Title", "New title:")
        if ok and t: self.title_text.set_text(t); self.draw()

    def edit_x_label(self):
        t, ok = QInputDialog.getText(self.parent(), "Edit X Label", "New x-axis:")
        if ok and t: self.x_label_text.set_text(t); self.draw()

    def edit_y_label(self):
        t, ok = QInputDialog.getText(self.parent(), "Edit Y Label", "New y-axis:")
        if ok and t: self.y_label_text.set_text(t); self.draw()

# ─── Main window ────────────────────────────────────────────────────────────

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
        self.current_vars = []

        # UI scaffolding
        main_w = QWidget(self); self.setCentralWidget(main_w)
        main_L = QHBoxLayout(main_w)

        # ─ Left panel ─────────────────────────────────────────────────────────
        left = QVBoxLayout()
        self.file_list = QListWidget()
        self.file_list.setStyleSheet("QListWidget { background:#2e2e2e; color:white; }")
        self.load_csv_files()
        self.file_list.itemClicked.connect(self.load_and_display_vars)
        left.addWidget(self.file_list)

        self.message_label = QLabel("")
        self.message_label.setStyleSheet("color:red;")
        left.addWidget(self.message_label)

        self.num_vars_label = QLabel("Number of vars (excl Time): 0")
        self.num_vars_label.setStyleSheet("color:white;")
        left.addWidget(self.num_vars_label)

        self.vars_container = QScrollArea()
        self.vars_container.setWidgetResizable(True)
        self.vars_widget = QWidget()
        self.vars_layout = QVBoxLayout(self.vars_widget)
        self.vars_container.setWidget(self.vars_widget)
        left.addWidget(self.vars_container)

        btn_un = QPushButton("Unselect All")
        btn_un.clicked.connect(self.unselect_all_vars)
        btn_un.setStyleSheet(
            "QPushButton { background:#3a3a3a; color:white; }"
            "QPushButton:hover { background:#505050; }"
        )
        left.addWidget(btn_un)

        self.slider_label = QLabel("Number of Time Steps: 1")
        self.slider_label.setStyleSheet("color:white;")
        left.addWidget(self.slider_label)

        self.time_slider = QSlider(Qt.Horizontal)
        self.time_slider.setMinimum(1)
        self.time_slider.valueChanged.connect(self.update_slider_label)
        self.time_slider.valueChanged.connect(self.plot_selected_vars)
        left.addWidget(self.time_slider)

        self.ignore_zero_checkbox = QCheckBox("Ignore zero values")
        self.ignore_zero_checkbox.setStyleSheet("color:white;")
        self.ignore_zero_checkbox.stateChanged.connect(self.load_and_display_vars)
        left.addWidget(self.ignore_zero_checkbox)

        self.ignore_dup_checkbox = QCheckBox("Ignore duplicates")
        self.ignore_dup_checkbox.setStyleSheet("color:white;")
        self.ignore_dup_checkbox.stateChanged.connect(self.load_and_display_vars)
        left.addWidget(self.ignore_dup_checkbox)

        self.default_name_checkbox = QCheckBox("Use default name")
        self.default_name_checkbox.setStyleSheet("color:white;")
        left.addWidget(self.default_name_checkbox)

        btn_save = QPushButton("Save Plot")
        btn_save.clicked.connect(self.save_plot)
        btn_save.setStyleSheet(
            "QPushButton { background:#3a3a3a; color:white; }"
            "QPushButton:hover { background:#505050; }"
        )
        left.addWidget(btn_save)
        # ─────────────────────────────────────────────────────────────────────

        # ─ Right panel ───────────────────────────────────────────────────────
        self.canvas = PlotCanvas(self)
        right = QVBoxLayout()
        right.addWidget(self.canvas)
        # ─────────────────────────────────────────────────────────────────────

        main_L.addLayout(left,1)
        main_L.addLayout(right,3)

    # ─── CSV loading ────────────────────────────────────────────────────────
    def load_csv_files(self):
        folder = os.path.join(os.getcwd(), "data", "csv")
        files  = glob.glob(os.path.join(folder, "*.csv"))
        files.sort(key=os.path.getmtime, reverse=True)
        for p in files:
            name = os.path.splitext(os.path.basename(p))[0]
            self.file_list.addItem(name)
            self.file_paths[name] = p

    def read_csv(self, name):
        with open(self.file_paths[name], "r") as f:
            reader = csv.reader(f)
            data   = [list(map(float,row)) for row in reader]
        time_list = data.pop(-1)
        return list(map(list, zip(*data))), time_list

    # ─── Populate checkboxes ────────────────────────────────────────────────
    def load_and_display_vars(self):
        itm = self.file_list.currentItem()
        if itm:
            self.current_file_name = itm.text()
            self.traced_flag = "traced" in itm.text().lower()
            self.current_spectrum, self.current_time_list = self.read_csv(self.current_file_name)
            n = len(self.current_time_list)
            self.time_slider.setMaximum(n)
            self.time_slider.setValue(n)
        elif not self.current_file_name:
            return

        # compute N from covariance rows
        total = len(self.current_spectrum)
        cov   = total - (3 if self.traced_flag else 0)
        raw   = cov / 48 if cov>0 else 0
        N     = int(round(raw**0.25)) if raw>0 else 0

        # clear old widgets & our list
        for i in reversed(range(self.vars_layout.count())):
            w = self.vars_layout.itemAt(i).widget()
            if w: w.setParent(None)
        self.series_checkboxes = []

        # “Time” checkbox (always shown)
        cb_t = QCheckBox("Time")
        cb_t.setStyleSheet("color:white;")
        cb_t.stateChanged.connect(self.plot_selected_vars)
        self.vars_layout.addWidget(cb_t)

        # now build one cb per series_index=i
        for i, series in enumerate(self.current_spectrum):
            if self.traced_flag and i<3:
                label = ["<<Tr(XX)>>/N","<<Tr(PP)>>/N","<<Tr(XP)>>/N"][i]
            else:
                cov_idx  = i - (3 if self.traced_flag else 0)
                num_blk  = 16 * (N**4)
                blk_idx  = cov_idx // num_blk
                blk_name = ['XX','PP','XP'][blk_idx]
                a1,i1,l1,m1,a2,i2,l2,m2 = reverse_indexXX(cov_idx, N)
                label = (f"<<{blk_name}>>: a1={a1}, i1={i1}, l1={l1}, m1={m1}, "
                         f"a2={a2}, i2={i2}, l2={l2}, m2={m2}")

            cb = QCheckBox(label)
            cb.setStyleSheet("color:white;")
            cb.series_index = i
            cb.stateChanged.connect(self.plot_selected_vars)
            self.vars_layout.addWidget(cb)
            self.series_checkboxes.append(cb)

        # now apply filtering by hiding unwanted boxes
        seen = set()
        for cb in self.series_checkboxes:
            idx  = cb.series_index
            data = self.current_spectrum[idx]
            hide = False

            # ignore all‐zero series?
            if self.ignore_zero_checkbox.isChecked() and all(v==0 for v in data):
                hide = True

            # ignore duplicates?
            if self.ignore_dup_checkbox.isChecked():
                key = tuple(round(v,4) for v in data)
                if key in seen:
                    hide = True
                else:
                    seen.add(key)

            cb.setVisible(not hide)

        # update count of visible (excluding Time)
        visible_count = sum(1 for cb in self.series_checkboxes if cb.isVisible())
        self.num_vars_label.setText(f"Number of vars (excl Time): {visible_count}")

    # ─── Helpers ────────────────────────────────────────────────────────────
    def unselect_all_vars(self):
        for cb in self.vars_widget.findChildren(QCheckBox):
            cb.setChecked(False)

    def update_slider_label(self, v):
        self.slider_label.setText(f"Number of Time Steps: {v}")
        self.selected_steps = v

    # ─── Build sel[], then plot ─────────────────────────────────────────────
    def plot_selected_vars(self):
        # compute N again
        total = len(self.current_spectrum)
        cov   = total - (3 if self.traced_flag else 0)
        raw   = cov / 48 if cov>0 else 0
        N     = int(round(raw**0.25)) if raw>0 else 0

        # legend map for traced
        traced_legend = {
            0: r'$\dfrac{<\mathrm{Tr}(X_i X_i)>}{N}$',
            1: r'$\dfrac{<\mathrm{Tr}(P_i P_i)>}{N}$',
            2: r'$\dfrac{<\mathrm{Tr}(X_i P_i)>}{N}$'
        }

        # collect selections by series_index
        sel = []
        for cb in self.vars_widget.findChildren(QCheckBox):
            if not cb.isChecked():
                continue
            if hasattr(cb, 'series_index'):
                sel.append(cb.series_index)
            else:
                # “Time”
                sel.append("Time")

        # time‐series path
        if "Time" in sel:
            if len(sel) < 2:
                self.message_label.setText("Select at least one variable in addition to Time.")
                return
            self.message_label.clear()

            x_axis  = self.current_time_list[:self.selected_steps]
            y_axes, y_titles = [], []

            for v in sel:
                if v == "Time": continue
                y_axes.append(self.current_spectrum[v][:self.selected_steps])
                if self.traced_flag and v < 3:
                    y_titles.append(traced_legend[v])
                else:
                    cov_idx = v - (3 if self.traced_flag else 0)
                    y_titles.append(generate_var_label(cov_idx, N))

            x_label = "Time"

        # scatter / pair‐plot path
        else:
            if len(sel) != 2:
                self.message_label.setText("Select exactly two variables for plotting.")
                return
            self.message_label.clear()

            i1, i2 = sel
            x_axis  = self.current_spectrum[i1][:self.selected_steps]
            y_axes  = [self.current_spectrum[i2][:self.selected_steps]]
            x_label = f"x{i1+1}"
            y_titles= [f"x{i2+1}"]

        # length check
        for y in y_axes:
            if len(y) != len(x_axis):
                self.message_label.setText("Error: x and y data have different lengths.")
                return

        title = f"Equations of Motion for {self.current_file_name}"
        self.canvas.plot(x_axis, y_axes, x_label, y_titles, title)
        self.current_vars = [x_label] + y_titles

    # ─── Save ────────────────────────────────────────────────────────────────
    def save_plot(self):
        if not self.current_vars:
            QMessageBox.warning(self, "Warning", "Please select variables first.")
            return
        first, *rest = self.current_vars
        if self.default_name_checkbox.isChecked():
            fname = f"{self.current_file_name}_{first}_{'_'.join(rest)}.png"
        else:
            nm, ok = QInputDialog.getText(self, "Save As", "Custom name:")
            if not ok or not nm: return
            fname = f"{nm}.png"
        out = os.path.join(os.getcwd(), "data", "plots", fname)
        self.canvas.figure.savefig(out, facecolor="#2e2e2e")
        QMessageBox.information(self, "Plot Saved", f"Saved: {out}")

# ─── run ───────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = MainWindow()
    w.show()
    sys.exit(app.exec_())
