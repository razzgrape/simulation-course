import tkinter as tk
from tkinter import ttk, messagebox
import math
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

class BallisticsApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Моделирование полета")
        self.root.geometry("1250x800")
        
        self.G = 9.81
        self.RHO = 1.225

        self.setup_ui()
        self.setup_plot()

    def setup_ui(self):
        control_frame = tk.Frame(self.root, width=420)
        control_frame.pack(side=tk.LEFT, fill=tk.Y, padx=10, pady=10)
        control_frame.pack_propagate(False)

        input_frame = tk.LabelFrame(control_frame, text=" Параметры снаряда ", padx=10, pady=10)
        input_frame.pack(side=tk.TOP, fill=tk.X, pady=5)

        params = [
            ("Нач. высота (м):", "0"),
            ("Угол (град):", "45"),
            ("Скорость (м/с):", "100"),
            ("Площадь (S):", "0.01"),
            ("Вес (кг):", "1.0"),
            ("Коэф. сопр. (C):", "0.47"),
            ("Шаг dt (с):", "0.01")
        ]
        
        self.entries = {}
        for i, (label, default) in enumerate(params):
            row, col = i // 2, (i % 2) * 2
            tk.Label(input_frame, text=label).grid(row=row, column=col, sticky="e", pady=2)
            entry = tk.Entry(input_frame, width=10)
            entry.insert(0, default)
            entry.grid(row=row, column=col+1, padx=5, pady=2)
            self.entries[label] = entry

        btn_frame = tk.Frame(control_frame)
        btn_frame.pack(side=tk.TOP, fill=tk.X, pady=10)
        
        tk.Button(btn_frame, text="РАССЧИТАТЬ", command=self.run_single, bg="#4CAF50", fg="white", width=18).grid(row=0, column=0, padx=5, pady=5)
        tk.Button(btn_frame, text="ЗАПУСТИТЬ ВСЕ ШАГИ", command=self.run_all_steps, bg="#2196F3", fg="white", width=18).grid(row=0, column=1, padx=5, pady=5)
        tk.Button(btn_frame, text="Очистить таблицу", command=self.clear_table, width=18).grid(row=1, column=0, padx=5, pady=5)
        tk.Button(btn_frame, text="Очистить график", command=self.clear_plot, width=18).grid(row=1, column=1, padx=5, pady=5)

        columns = ("step", "dist", "height", "speed")
        self.tree = ttk.Treeview(control_frame, columns=columns, show="headings")
        self.tree.heading("step", text="Шаг, с")
        self.tree.heading("dist", text="Дальность, м")
        self.tree.heading("height", text="Высота, м")
        self.tree.heading("speed", text="V кон, м/с")
        self.tree.pack(side=tk.TOP, fill=tk.BOTH, expand=True, pady=10)

        for col in columns:
            self.tree.column(col, width=85, anchor=tk.CENTER)

    def setup_plot(self):
        plot_frame = tk.Frame(self.root)
        plot_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=10, pady=10)
        self.fig, self.ax = plt.subplots()
        
        self.ax.set_title("Траектории полета тела в атмосфере")
        self.ax.set_xlabel("Дальность полета (метры)")
        self.ax.set_ylabel("Высота полета (метры)")
        self.ax.grid(True)
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    def simulate(self, dt_val):
        y0 = float(self.entries["Нач. высота (м):"].get())
        alpha = math.radians(float(self.entries["Угол (град):"].get()))
        v0 = float(self.entries["Скорость (м/с):"].get())
        s_val = float(self.entries["Площадь (S):"].get())
        m = float(self.entries["Вес (кг):"].get())
        c_val = float(self.entries["Коэф. сопр. (C):"].get())

        k = (c_val * s_val * self.RHO) / (2 * m)

        x, y = 0.0, y0
        vx = v0 * math.cos(alpha)
        vy = v0 * math.sin(alpha)
        max_h = y0
        
        x_pts, y_pts = [x], [y]

        while y > 0 or (len(x_pts) == 1): 
            v = math.sqrt(vx**2 + vy**2)
            if y > max_h: max_h = y
            
            vx_next = vx - k * vx * v * dt_val
            vy_next = vy - (self.G + k * vy * v) * dt_val
            
            y_next = y + vy_next * dt_val
            
            if y_next < 0:
                # Если падает, вычисляем долю шага до y=0
                fraction = y / (y - y_next)
                x += vx_next * dt_val * fraction
                y = 0.0
                vx, vy = vx_next, vy_next
                x_pts.append(x)
                y_pts.append(y)
                break
            
            x += vx_next * dt_val
            y = y_next
            vx, vy = vx_next, vy_next

            x_pts.append(x)
            y_pts.append(y)
            
        return x_pts, y_pts, max_h, math.sqrt(vx**2 + vy**2)

    def run_single(self):
        try:
            dt = float(self.entries["Шаг dt (с):"].get())
            self.execute_and_plot(dt)
        except Exception as e:
            messagebox.showerror("Ошибка", f"Проверьте ввод: {e}")

    def run_all_steps(self):
        steps = [1.0, 0.1, 0.01, 0.001, 0.0001]
        for dt in steps:
            self.execute_and_plot(dt)

    def execute_and_plot(self, dt):
        x_pts, y_pts, max_h, v_final = self.simulate(dt)
        self.ax.plot(x_pts, y_pts, label=f"Шаг dt={dt}с")
        self.ax.legend()
        self.ax.relim()
        self.ax.autoscale_view()
        self.canvas.draw()
        self.tree.insert("", tk.END, values=(f"{dt}", f"{x_pts[-1]:.2f}", f"{max_h:.2f}", f"{v_final:.2f}"))
        self.root.update()

    def clear_table(self):
        for item in self.tree.get_children(): self.tree.delete(item)

    def clear_plot(self):
        self.ax.clear()
        self.ax.set_title("Траектории полета тела в атмосфере")
        self.ax.set_xlabel("Дальность полета (метры)")
        self.ax.set_ylabel("Высота полета (метры)")
        self.ax.grid(True)
        self.canvas.draw()

if __name__ == "__main__":
    root = tk.Tk()
    app = BallisticsApp(root)
    root.mainloop()