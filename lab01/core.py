import sys
import math
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                             QHBoxLayout, QFormLayout, QLineEdit, QPushButton, 
                             QTreeWidget, QTreeWidgetItem, QGroupBox, QMessageBox)

class BallisticsApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Баллистический симулятор")
        self.setMinimumSize(1200, 800)
        
        self.G = 9.81
        self.RHO = 1.225
        self.color_map = {
            1.0: '#007AFF',    
            0.1: '#FF9500',    
            0.01: '#34C759',   
            0.001: '#FF3B30',  
            0.0001: '#000000'  
        }

        self.setup_ui()

    def setup_ui(self):
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QHBoxLayout(central_widget)

        side_panel = QVBoxLayout()
        side_panel.setContentsMargins(10, 10, 10, 10)
        
        params_group = QGroupBox("Параметры снаряда")
        form_layout = QFormLayout()
        
        self.inputs = {}
        fields = [
            ("Нач. высота (м):", "0"), ("Угол (град):", "45"),
            ("Скорость (м/с):", "100"), ("Площадь (S):", "0.01"),
            ("Вес (кг):", "1.0"), ("Коэф. сопр. (C):", "0.47"),
            ("Шаг dt (с):", "0.01")
        ]
        
        for label_text, default in fields:
            input_field = QLineEdit(default)
            form_layout.addRow(label_text, input_field)
            self.inputs[label_text] = input_field
            
        params_group.setLayout(form_layout)
        side_panel.addWidget(params_group)

        self.btn_calc = QPushButton("РАССЧИТАТЬ")
        self.btn_calc.clicked.connect(self.run_single)
        self.btn_all = QPushButton("ЗАПУСТИТЬ ВСЕ ШАГИ")
        self.btn_all.clicked.connect(self.run_all_steps)
        self.btn_clear = QPushButton("ОЧИСТИТЬ ВСЕ")
        self.btn_clear.clicked.connect(self.clear_all)

        side_panel.addWidget(self.btn_calc)
        side_panel.addWidget(self.btn_all)
        side_panel.addWidget(self.btn_clear)

        self.tree = QTreeWidget()
        self.tree.setHeaderLabels(["Шаг, с", "Дальность", "Высота", "V кон"])
        self.tree.setColumnWidth(0, 60)
        side_panel.addWidget(self.tree)

        main_layout.addLayout(side_panel, 1)

        self.figure, self.ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
        self.canvas = FigureCanvas(self.figure)
        
        self.ax.set_title("Траектория полета", fontsize=14, fontweight='bold')
        self.ax.set_xlabel("Расстояние (м)")
        self.ax.set_ylabel("Высота (м)")
        self.ax.grid(True, linestyle=':', alpha=0.6)
        
        main_layout.addWidget(self.canvas, 3)

    def get_float_val(self, key):
        text = self.inputs[key].text().replace(',', '.').strip()
        return float(text)

    def simulate(self, dt_val):
        try:
            y0 = self.get_float_val("Нач. высота (м):")
            alpha = math.radians(self.get_float_val("Угол (град):"))
            v0 = self.get_float_val("Скорость (м/с):")
            s_val = self.get_float_val("Площадь (S):")
            m = self.get_float_val("Вес (кг):")
            c_val = self.get_float_val("Коэф. сопр. (C):")
        except ValueError:
            QMessageBox.warning(self, "Ошибка ввода", "Пожалуйста, введите корректные числа во все поля!")
            return None

        k = (c_val * s_val * self.RHO) / (2 * m)
        x, y = 0.0, y0
        vx, vy = v0 * math.cos(alpha), v0 * math.sin(alpha)
        max_h = y
        x_pts, y_pts = [x], [y]

        max_steps = 500000 
        steps = 0

        while (y > 0 or len(x_pts) == 1) and steps < max_steps:
            steps += 1
            v = math.sqrt(vx**2 + vy**2)
            if y > max_h: max_h = y
            
            vx_next = vx - k * vx * v * dt_val
            vy_next = vy - (self.G + k * vy * v) * dt_val
            y_next = y + vy_next * dt_val
            
            if y_next < 0:
                fraction = y / (y - y_next) if (y - y_next) != 0 else 0
                x += vx_next * dt_val * fraction
                y = 0.0
                x_pts.append(x); y_pts.append(y)
                break
            
            x += vx_next * dt_val
            y = y_next
            vx, vy = vx_next, vy_next
            x_pts.append(x); y_pts.append(y)
            
        if steps >= max_steps:
             QMessageBox.warning(self, "Сбой физики", "Снаряд улетел слишком далеко. Проверьте параметры!")
             return None

        return x_pts, y_pts, max_h, math.sqrt(vx**2 + vy**2)

    def execute_and_plot(self, dt):
        res = self.simulate(dt)
        if not res: 
            return 
        
        x_pts, y_pts, max_h, v_final = res
        color = self.color_map.get(dt, '#5856D6') 
        
        self.ax.plot(x_pts, y_pts, label=f"dt={dt}s", color=color, linewidth=2)
        
        self.ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1), frameon=False, fontsize=10)
        
        self.ax.relim()
        self.ax.autoscale_view()
        self.canvas.draw()
        
        item = QTreeWidgetItem([str(dt), f"{x_pts[-1]:.2f}", f"{max_h:.2f}", f"{v_final:.2f}"])
        self.tree.addTopLevelItem(item)

    def run_single(self):
        try:
            dt = self.get_float_val("Шаг dt (с):")
            self.execute_and_plot(dt)
        except ValueError:
            QMessageBox.warning(self, "Ошибка", "Проверьте поле 'Шаг dt (с)'. Там должно быть число.")

    def run_all_steps(self):
        self.clear_all() 
        for dt in [1.0, 0.1, 0.01, 0.001, 0.0001]:
            self.execute_and_plot(dt)

    def clear_all(self):
        self.ax.clear()
        self.ax.set_title("Траектория полета", fontsize=14, fontweight='bold')
        self.ax.set_xlabel("Расстояние (м)")
        self.ax.set_ylabel("Высота (м)")
        self.ax.grid(True, linestyle=':', alpha=0.6)
        self.canvas.draw()
        self.tree.clear()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    app.setStyle("Fusion") 
    window = BallisticsApp()
    window.show()
    sys.exit(app.exec())