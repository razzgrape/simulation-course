import customtkinter as ctk
import numpy as np
import math
from collections import Counter
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

from core import simulate_poisson_stream, get_theoretical_probs, calculate_empirical_stats

ctk.set_appearance_mode("dark")
ctk.set_default_color_theme("blue")

ACCENT  = "#4F8EF7"
ACCENT2 = "#F76B4F"
BG      = "#0F1117"
CARD    = "#1A1D27"
CARD2   = "#21253A"
TEXT    = "#E8EAF0"
MUTED   = "#6B7280"
SUCCESS = "#34D399"


class StatCard(ctk.CTkFrame):
    def __init__(self, master, label, **kwargs):
        super().__init__(master, fg_color=CARD2, corner_radius=14, **kwargs)
        ctk.CTkLabel(self, text=label, font=("Helvetica", 11), text_color=MUTED).pack(pady=(14, 2))
        self.value_lbl = ctk.CTkLabel(self, text="—", font=("Helvetica", 26, "bold"), text_color=TEXT)
        self.value_lbl.pack(pady=(0, 14))

    def set(self, val):
        self.value_lbl.configure(text=val)


class ParamRow(ctk.CTkFrame):
    def __init__(self, master, label, default, **kwargs):
        super().__init__(master, fg_color="transparent", **kwargs)
        
        # Заголовок параметра
        ctk.CTkLabel(self, text=label, font=("Helvetica", 12), text_color=MUTED,
                     width=200, anchor="w").grid(row=0, column=0, padx=(0, 12), sticky="w")

        # Поле ввода вместо слайдера
        self.entry = ctk.CTkEntry(self, width=120, fg_color=CARD2, 
                                  border_color="#2D3147", text_color=ACCENT,
                                  font=("Helvetica", 13, "bold"))
        self.entry.insert(0, str(default))  # Установка значения по умолчанию
        self.entry.grid(row=0, column=1, padx=(0, 14))

    def get(self):
        try:
            # Превращаем текст из поля в число
            return int(self.entry.get())
        except ValueError:
            # Если введено не число, возвращаем 0 или можно выводить ошибку
            return 0


class App(ctk.CTk):
    def __init__(self):
        super().__init__()
        self.title("Пуассоновский поток · Server Simulation")
        self.geometry("980x760")
        self.configure(fg_color=BG)
        self.resizable(True, True)
        self._build()

    def _build(self):
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(2, weight=1)
        self._build_header()
        self._build_controls()
        self._build_main()
        self._build_conclusion()

    def _build_header(self):
        hdr = ctk.CTkFrame(self, fg_color=CARD, corner_radius=0, height=64)
        hdr.grid(row=0, column=0, sticky="ew")
        hdr.grid_propagate(False)
        hdr.grid_columnconfigure(0, weight=1)
        ctk.CTkLabel(hdr, text="  ◈  Poisson Server Flow Simulator",
                     font=("Helvetica", 18, "bold"), text_color=TEXT,
                     anchor="w").grid(row=0, column=0, padx=24, pady=16, sticky="w")
        ctk.CTkLabel(hdr, text="Моделирование простейшего потока событий",
                     font=("Helvetica", 11), text_color=MUTED,
                     anchor="e").grid(row=0, column=1, padx=24, pady=16, sticky="e")

    def _build_controls(self):
        panel = ctk.CTkFrame(self, fg_color=CARD, corner_radius=18)
        panel.grid(row=1, column=0, padx=20, pady=(18, 10), sticky="ew")
        panel.grid_columnconfigure(0, weight=1)

        inner = ctk.CTkFrame(panel, fg_color="transparent")
        inner.grid(row=0, column=0, padx=28, pady=18, sticky="ew")

        # Теперь передаем только название и значение по умолчанию
        self.p_lambda = ParamRow(inner, "Интенсивность  λ  (зап/с)", 5)
        self.p_lambda.grid(row=0, column=0, pady=6, sticky="w")

        self.p_T = ParamRow(inner, "Интервал  T  (секунд)", 3)
        self.p_T.grid(row=1, column=0, pady=6, sticky="w")

        self.p_N = ParamRow(inner, "Число прогонов  N", 1000)
        self.p_N.grid(row=2, column=0, pady=6, sticky="w")

        self.btn = ctk.CTkButton(inner, text="▶   Запустить симуляцию",
                                  font=("Helvetica", 13, "bold"),
                                  fg_color=ACCENT, hover_color="#3A72D8",
                                  text_color="white", corner_radius=12,
                                  height=42, width=240, command=self.run)
        self.btn.grid(row=0, column=1, rowspan=3, padx=(40, 0))

    def _build_main(self):
        main = ctk.CTkFrame(self, fg_color="transparent")
        main.grid(row=2, column=0, padx=20, pady=0, sticky="nsew")
        main.grid_columnconfigure(0, weight=3)
        main.grid_columnconfigure(1, weight=1)
        main.grid_rowconfigure(0, weight=1)
        self._build_chart(main)
        self._build_stats(main)

    def _build_chart(self, parent):
        chart_frame = ctk.CTkFrame(parent, fg_color=CARD, corner_radius=18)
        chart_frame.grid(row=0, column=0, padx=(0, 10), pady=0, sticky="nsew")
        chart_frame.grid_rowconfigure(0, weight=1)
        chart_frame.grid_columnconfigure(0, weight=1)

        self.fig = Figure(facecolor=CARD)
        self.ax = self.fig.add_subplot(111, facecolor=CARD)
        self._style_ax_empty()

        self.canvas = FigureCanvasTkAgg(self.fig, master=chart_frame)
        self.canvas.get_tk_widget().configure(bg=CARD)
        self.canvas.get_tk_widget().grid(row=0, column=0, padx=18, pady=18, sticky="nsew")

    def _style_ax_empty(self):
        ax = self.ax
        ax.set_title("Запустите симуляцию", color=MUTED, fontsize=12, pad=14)
        ax.set_xlabel("Число запросов k", color=MUTED, fontsize=10)
        ax.set_ylabel("Вероятность  P(X = k)", color=MUTED, fontsize=10)
        ax.tick_params(colors=MUTED)
        for spine in ax.spines.values():
            spine.set_color("#2D3147")
            spine.set_linewidth(0.8)
        ax.grid(axis="y", color="#2D3147", linewidth=0.8, linestyle="--")

    def _build_stats(self, parent):
        col = ctk.CTkFrame(parent, fg_color="transparent")
        col.grid(row=0, column=1, sticky="nsew")
        col.grid_columnconfigure(0, weight=1)

        for i, (lbl, attr) in enumerate([
            ("μ = λT  (теор.)",    "card_mu"),
            ("Среднее  (эмп.)",    "card_mean"),
            ("Дисперсия  (эмп.)",  "card_var"),
            ("Дисперсия  (теор.)", "card_varth"),
        ]):
            c = StatCard(col, lbl)
            c.grid(row=i, column=0, pady=6, sticky="ew")
            setattr(self, attr, c)

    def _build_conclusion(self):
        box = ctk.CTkFrame(self, fg_color=CARD, corner_radius=18)
        box.grid(row=3, column=0, padx=20, pady=(10, 20), sticky="ew")
        box.grid_columnconfigure(1, weight=1)
        ctk.CTkLabel(box, text="◆", font=("Helvetica", 14),
                     text_color=ACCENT).grid(row=0, column=0, padx=(18, 8), pady=16)
        self.conclusion = ctk.CTkLabel(
            box, text="После запуска симуляции здесь появится вывод.",
            font=("Helvetica", 12), text_color=MUTED,
            anchor="w", wraplength=860, justify="left")
        self.conclusion.grid(row=0, column=1, padx=(0, 18), pady=16, sticky="ew")

    def run(self):
        self.btn.configure(state="disabled", text="⏳  Считаем...")
        self.update()

        lmbda = self.p_lambda.get()
        T     = self.p_T.get()
        N     = self.p_N.get()
        mu    = lmbda * T

        all_counts = [len(simulate_poisson_stream(lmbda, T)) for _ in range(N)]
        emp_mean, emp_var, freqs = calculate_empirical_stats(all_counts)

        self.card_mu.set(f"{mu:.2f}")
        self.card_mean.set(f"{emp_mean:.2f}")
        self.card_var.set(f"{emp_var:.2f}")
        self.card_varth.set(f"{mu:.2f}")

        self._draw(freqs, N, mu, lmbda, T)
        self._conclude(lmbda, T, N, mu, emp_mean, emp_var)
        self.btn.configure(state="normal", text="▶   Запустить симуляцию")

    def _draw(self, freqs, N, mu, lmbda, T):
        self.ax.clear()

        max_k = max(freqs.keys())
        ks = list(range(max_k + 2))
        emp = [freqs.get(k, 0) / N for k in ks]
        th_dict = get_theoretical_probs(mu, 1, max_k + 1)
        th = [th_dict.get(k, 0.0) for k in ks]

        self.ax.bar(ks, emp, width=0.55, color=ACCENT, alpha=0.72,
                    label="Эмпирическое", zorder=2)
        self.ax.plot(ks, th, color=ACCENT2, linewidth=2.2, marker="o",
                     markersize=5, label="Пуассон (теор.)", zorder=3)

        self.ax.set_title(
            f"Распределение запросов   λ={lmbda}  T={T}с  N={N}",
            color=TEXT, fontsize=12, pad=14)
        self.ax.set_xlabel("Число запросов k", color=MUTED, fontsize=10)
        self.ax.set_ylabel("Вероятность  P(X = k)", color=MUTED, fontsize=10)
        self.ax.tick_params(colors=MUTED)
        for spine in self.ax.spines.values():
            spine.set_color("#2D3147")
            spine.set_linewidth(0.8)
        self.ax.grid(axis="y", color="#2D3147", linewidth=0.8, linestyle="--")
        self.ax.legend(fontsize=10, facecolor=CARD2, edgecolor="#2D3147",
                       labelcolor=TEXT, framealpha=1)
        self.canvas.draw()

    def _conclude(self, lmbda, T, N, mu, emp_mean, emp_var):
        diff_m = abs(emp_mean - mu)
        diff_v = abs(emp_var - mu)
        ok = diff_m < 0.5 and diff_v < 1.0
        verdict = "Результаты хорошо согласуются с теорией ✓" if ok else \
                  "Увеличьте N для лучшего совпадения с теорией."
        text = (
            f"Смоделировано {N} интервалов T={T}с при λ={lmbda} зап/с.  "
            f"Теоретическое μ = λT = {mu:.2f}.  "
            f"Эмп. среднее = {emp_mean:.2f} (Δ {diff_m:.2f}),  "
            f"эмп. дисперсия = {emp_var:.2f} (Δ от теор. {diff_v:.2f}).  "
            f"Для распр. Пуассона M[X] = D[X] = μ.  {verdict}"
        )
        self.conclusion.configure(text=text, text_color=SUCCESS if ok else MUTED)


if __name__ == "__main__":
    app = App()
    app.mainloop()