import tkinter as tk
from tkinter import messagebox
import numpy as np
import time, threading
import matplotlib

matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.ticker import AutoMinorLocator
import matplotlib.patheffects as pe

from solver import solve_heat


# ─── Цвета ───────────────────────────────────────────────────────────────────

BG           = "#0d1117"
PANEL        = "#161b22"
PANEL2       = "#1c2333"
BORDER       = "#30363d"
ACCENT       = "#58a6ff"
ACCENT_HOT   = "#ff7b72"
ACCENT_GOLD  = "#e3b341"
TEXT         = "#e6edf3"
TEXT_DIM     = "#8b949e"
GREEN        = "#3fb950"
RED          = "#f85149"

TABLE_TAUS   = [0.1, 0.01, 0.001, 0.0001]
TABLE_HS     = [0.1, 0.01, 0.001, 0.0001]

class HeatApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Теплопроводность · МКР · Метод прогонки")
        self.configure(bg=BG)
        self.geometry("1490x920")
        self.minsize(1100, 720)
        self._build()

    def _build(self):
        self.columnconfigure(0, weight=0, minsize=330)
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=1)
        self._panel_left()
        self._panel_right()

    def _panel_left(self):
        lf = tk.Frame(self, bg=PANEL, width=330)
        lf.grid(row=0, column=0, sticky="nsew", padx=(14, 0), pady=14)
        lf.grid_propagate(False)
        lf.columnconfigure(1, weight=1)

        r = 0
        hdr = tk.Frame(lf, bg=PANEL)
        hdr.grid(row=r, column=0, columnspan=2, sticky="ew",
                 padx=16, pady=(20, 6))
        tk.Label(hdr, text="ρc", bg=PANEL, fg=TEXT_DIM,
                 font=("Georgia", 14, "italic")).pack(side="left", pady=2)
        tk.Label(hdr, text=" ∂T/∂t ", bg=PANEL, fg=ACCENT,
                 font=("Georgia", 18, "bold italic")).pack(side="left")
        tk.Label(hdr, text="= λ ∂²T/∂x²", bg=PANEL, fg=TEXT_DIM,
                 font=("Georgia", 14, "italic")).pack(side="left", pady=2)
        r += 1

        tk.Frame(lf, bg=BORDER, height=1).grid(
            row=r, column=0, columnspan=2, sticky="ew", padx=10, pady=(0, 10))
        r += 1

        r = self._section_row(lf, r, "МАТЕРИАЛ")
        for label, key, dflt in [
            ("ρ,  кг/м³",       "rho",    "8960"),
            ("c,  Дж/(кг·К)",  "c_heat", "385"),
            ("λ,  Вт/(м·К)",   "lam",    "385"),
            ("L,  м",          "L",      "0.05"),
        ]:
            r = self._row(lf, r, label, key, dflt)

        r = self._section_row(lf, r, "УСЛОВИЯ")
        for label, key, dflt in [
            ("T₀,  °C",        "T0",      "20"),
            ("T_лев,  °C",     "T_left",  "200"),
            ("T_пр,  °C",      "T_right", "50"),
        ]:
            r = self._row(lf, r, label, key, dflt)

        r = self._section_row(lf, r, "ШАГИ")
        for label, key, dflt in [
            ("τ,  с",          "tau",    "0.001"),
            ("h,  м",          "h",      "0.001"),
            ("t_end,  с",      "t_end",  "2.0"),
        ]:
            r = self._row(lf, r, label, key, dflt)

        tk.Frame(lf, bg=BORDER, height=1).grid(
            row=r, column=0, columnspan=2, sticky="ew", padx=10, pady=10)
        r += 1

        self.btn1 = self._mkbtn(lf, r, "▶  Запуск",          ACCENT,      self._run_single); r += 1
        self.btn2 = self._mkbtn(lf, r, "⏩  Заполнить таблицу", ACCENT_HOT, self._run_table);  r += 1

        tk.Frame(lf, bg=BORDER, height=1).grid(
            row=r, column=0, columnspan=2, sticky="ew", padx=10, pady=8)
        r += 1

        for txt, var_attr, color, font_spec in [
            ("Температура в центре:", "var_tc",   GREEN,    ("Courier New", 16, "bold")),
            ("Время расчёта:",        "var_time", TEXT,     ("Courier New", 10)),
        ]:
            tk.Label(lf, text=txt, bg=PANEL, fg=TEXT_DIM,
                     font=("Courier New", 8)).grid(
                row=r, column=0, columnspan=2, sticky="w", padx=18)
            r += 1
            var = tk.StringVar(value="—")
            setattr(self, var_attr, var)
            tk.Label(lf, textvariable=var, bg=PANEL, fg=color,
                     font=font_spec).grid(
                row=r, column=0, columnspan=2, sticky="w", padx=18)
            r += 1

        self.var_status = tk.StringVar(value="Готов к расчёту")
        tk.Label(lf, textvariable=self.var_status, bg=PANEL, fg=TEXT_DIM,
                 font=("Courier New", 8), wraplength=300,
                 justify="left").grid(
            row=r, column=0, columnspan=2, sticky="w", padx=18, pady=(8, 6))

    def _section_row(self, parent, row, text):
        f = tk.Frame(parent, bg=PANEL2)
        f.grid(row=row, column=0, columnspan=2, sticky="ew",
               padx=10, pady=(8, 2))
        tk.Label(f, text=f"  {text}", bg=PANEL2, fg=ACCENT_GOLD,
                 font=("Courier New", 8, "bold")).pack(anchor="w", ipady=3)
        return row + 1

    def _row(self, parent, row, label, key, default):
        tk.Label(parent, text=label, bg=PANEL, fg=TEXT,
                 font=("Courier New", 10), anchor="w").grid(
            row=row, column=0, sticky="w", padx=(18, 4), pady=2)
        e = tk.Entry(parent, font=("Courier New", 11),
                     bg=PANEL2, fg=ACCENT, insertbackground=ACCENT,
                     relief="flat", highlightthickness=1,
                     highlightcolor=ACCENT, highlightbackground=BORDER,
                     width=11)
        e.insert(0, default)
        e.grid(row=row, column=1, sticky="ew", padx=(0, 14), pady=2)
        self.entries[key] = e
        return row + 1

    def _mkbtn(self, parent, row, text, color, cmd):
        b = tk.Button(parent, text=text,
                      font=("Courier New", 11, "bold"),
                      bg=color, fg=BG,
                      activebackground=color, activeforeground=BG,
                      relief="flat", cursor="hand2", command=cmd)
        b.grid(row=row, column=0, columnspan=2, sticky="ew",
               padx=12, ipady=9, pady=3)
        return b

    entries: dict = {}

    def _panel_right(self):
        rf = tk.Frame(self, bg=BG)
        rf.grid(row=0, column=1, sticky="nsew", padx=14, pady=14)
        rf.columnconfigure(0, weight=1)
        rf.rowconfigure(0, weight=1)
        rf.rowconfigure(1, weight=0)

        plot_frame = tk.Frame(rf, bg=PANEL,
                              highlightthickness=1,
                              highlightbackground=BORDER)
        plot_frame.grid(row=0, column=0, sticky="nsew", pady=(0, 10))

        self.fig = Figure(facecolor=PANEL)
        self.fig.subplots_adjust(left=0.07, right=0.97,
                                 top=0.89, bottom=0.11)
        self.ax = self.fig.add_subplot(111, facecolor="#060d14")
        self._style_ax()

        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.get_tk_widget().configure(bg=PANEL, highlightthickness=0)
        self.canvas.get_tk_widget().pack(fill="both", expand=True,
                                         padx=2, pady=2)
        self.canvas.draw()

        tbl = tk.Frame(rf, bg=PANEL,
                       highlightthickness=1,
                       highlightbackground=BORDER)
        tbl.grid(row=1, column=0, sticky="ew")

        hdr_f = tk.Frame(tbl, bg=PANEL2)
        hdr_f.pack(fill="x")
        tk.Label(hdr_f,
                 text="  Таблица сходимости:  T_центр (°C)  при  t = 2 с",
                 bg=PANEL2, fg=ACCENT_GOLD,
                 font=("Courier New", 10, "bold")).pack(
            side="left", pady=7, padx=8)
        tk.Label(hdr_f,
                 text="🟢 < 0.5%  🟡 < 2%  🔴 > 2%  (отклонение от наилучшего)",
                 bg=PANEL2, fg=TEXT_DIM,
                 font=("Courier New", 8)).pack(side="right", padx=12)

        grid_f = tk.Frame(tbl, bg=BG)
        grid_f.pack(fill="x", padx=6, pady=(2, 8))

        SH = dict(bg=PANEL2, fg=ACCENT,
                  font=("Courier New", 9, "bold"),
                  relief="flat", padx=14, pady=6)
        SV = dict(bg=PANEL, fg=TEXT_DIM,
                  font=("Courier New", 10),
                  relief="flat", padx=14, pady=5, width=13)

        tk.Label(grid_f, text="τ  \\  h", **SH).grid(
            row=0, column=0, padx=1, pady=1, sticky="ew")
        for ci, hv in enumerate(TABLE_HS):
            tk.Label(grid_f, text=f"h = {hv}", **SH).grid(
                row=0, column=ci + 1, padx=1, pady=1, sticky="ew")

        self._cvars = {}
        self._clbls = {}
        for ri, tv in enumerate(TABLE_TAUS):
            tk.Label(grid_f, text=f"τ = {tv}", **SH).grid(
                row=ri + 1, column=0, padx=1, pady=1, sticky="ew")
            for ci, hv in enumerate(TABLE_HS):
                var = tk.StringVar(value="")
                lbl = tk.Label(grid_f, textvariable=var, **SV)
                lbl.grid(row=ri + 1, column=ci + 1,
                         padx=1, pady=1, sticky="ew")
                self._cvars[(tv, hv)] = var
                self._clbls[(tv, hv)] = lbl

    def _style_ax(self):
        ax = self.ax
        ax.set_facecolor("#060d14")
        for sp in ax.spines.values():
            sp.set_color(BORDER)
        ax.tick_params(colors=TEXT_DIM, labelsize=9,
                       which="both", direction="in", top=True, right=True)
        ax.xaxis.set_minor_locator(AutoMinorLocator(5))
        ax.yaxis.set_minor_locator(AutoMinorLocator(5))
        ax.grid(which="major", color="#1a2535", linewidth=0.7)
        ax.grid(which="minor", color="#111c28", linewidth=0.4)
        ax.set_xlabel("x,  м",  color=TEXT_DIM, fontsize=10,
                      fontfamily="Courier New")
        ax.set_ylabel("T,  °C", color=TEXT_DIM, fontsize=10,
                      fontfamily="Courier New")
        ax.set_title("Распределение температуры вдоль пластины",
                     color=TEXT, fontsize=11, fontfamily="Courier New",
                     pad=10, loc="left")

    def _get_params(self):
        p = {}
        for k, e in self.entries.items():
            raw = e.get().strip()
            try:
                p[k] = float(raw)
            except ValueError:
                messagebox.showerror("Ошибка ввода",
                                     f"Неверное значение «{k}»: '{raw}'")
                return None
        return p

    def _lock(self):
        self.btn1.config(state="disabled")
        self.btn2.config(state="disabled")

    def _unlock(self):
        self.btn1.config(state="normal")
        self.btn2.config(state="normal")

    def _run_single(self):
        p = self._get_params()
        if p is None:
            return
        self._lock()
        self.var_status.set("Расчёт…")
        threading.Thread(target=self._th_single, args=(p,),
                         daemon=True).start()

    def _th_single(self, p):
        t0 = time.perf_counter()
        try:
            T, h_e, N = solve_heat(
                L=p["L"], T0=p["T0"],
                T_left=p["T_left"], T_right=p["T_right"],
                rho=p["rho"], c=p["c_heat"], lam=p["lam"],
                tau=p["tau"], h=p["h"], t_end=p["t_end"])
            dt = time.perf_counter() - t0
            self.after(0, self._done_single, T, h_e, N, dt, p)
        except Exception as ex:
            self.after(0, lambda: messagebox.showerror("Ошибка", str(ex)))
            self.after(0, lambda: self.var_status.set("Ошибка ✗"))
            self.after(0, self._unlock)

    def _done_single(self, T, h_e, N, dt, p):
        mid = N // 2
        Tc  = T[mid]
        self.var_tc.set(f"{Tc:.5f} °C")
        self.var_time.set(f"{dt * 1000:.2f} мс")
        self.var_status.set(
            f"N = {N}  ·  шагов по t = {round(p['t_end']/p['tau'])}  ·  h_эфф = {h_e:.4g} м")

        x  = np.linspace(0, p["L"], N)
        ax = self.ax
        ax.cla()
        self._style_ax()

        T_base = p["T0"]
        ax.fill_between(x, T, T_base,
                        where=(T >= T_base),
                        alpha=0.15, color=ACCENT_HOT, linewidth=0)
        ax.fill_between(x, T, T_base,
                        where=(T < T_base),
                        alpha=0.10, color=ACCENT,     linewidth=0)

        ax.plot(x, T, color=ACCENT_HOT, linewidth=2.8, zorder=4,
                path_effects=[
                    pe.SimpleLineShadow(offset=(0, -1.5),
                                        shadow_color=ACCENT_HOT,
                                        alpha=0.35, rho=1.0),
                    pe.Normal()
                ])

        ax.scatter([x[0], x[-1]], [T[0], T[-1]],
                   color=ACCENT_GOLD, s=70, zorder=7)

        ax.axvline(x[mid], color=TEXT_DIM, lw=0.9,
                   linestyle=":", alpha=0.5, zorder=3)
        ax.scatter([x[mid]], [Tc], color=GREEN, zorder=8, s=90)

        x_off = p["L"] * 0.04
        ax.annotate(
            f"  x = {x[mid]:.4g} м\n  T = {Tc:.5g} °C",
            xy=(x[mid], Tc),
            xytext=(x[mid] + x_off, Tc),
            color=GREEN, fontsize=9, fontfamily="Courier New",
            va="center",
            arrowprops=dict(arrowstyle="-", color=GREEN,
                            lw=0.8, alpha=0.5)
        )

        ax.text(0.01, 0.97,
                f"τ = {p['tau']:.4g} с  |  h = {p['h']:.4g} м  |  "
                f"t_end = {p['t_end']:.4g} с  |  N = {N}",
                transform=ax.transAxes,
                color=TEXT_DIM, fontsize=8,
                fontfamily="Courier New", va="top")

        self.canvas.draw()
        self._unlock()

    def _run_table(self):
        p = self._get_params()
        if p is None:
            return
        self._lock()
        for var in self._cvars.values():
            var.set("…")
        for lbl in self._clbls.values():
            lbl.configure(fg=TEXT_DIM)
        self.var_status.set("Заполняю таблицу…")
        threading.Thread(target=self._th_table, args=(p,),
                         daemon=True).start()

    def _th_table(self, p):
        total   = len(TABLE_TAUS) * len(TABLE_HS)
        done    = 0
        results = {}

        for tv in TABLE_TAUS:
            for hv in TABLE_HS:
                try:
                    T, _, N = solve_heat(
                        L=p["L"], T0=p["T0"],
                        T_left=p["T_left"], T_right=p["T_right"],
                        rho=p["rho"], c=p["c_heat"], lam=p["lam"],
                        tau=tv, h=hv, t_end=2.0)
                    val = T[N // 2]
                    results[(tv, hv)] = val
                    val_str = f"{val:.4f}"
                    err = False
                except Exception:
                    val_str = "ERR"
                    err = True

                done += 1

                def _upd(tv=tv, hv=hv, vs=val_str, e=err, d=done):
                    self._cvars[(tv, hv)].set(vs)
                    self._clbls[(tv, hv)].configure(
                        fg=RED if e else TEXT_DIM)
                    self.var_status.set(
                        f"Таблица: {d}/{total}  (τ={tv}, h={hv})")

                self.after(0, _upd)

        self.after(0, self._color_converged, results)
        self.after(0, lambda: self.var_status.set(
            f"Готово — {total} ячеек ✓  |  "
            "🟢 сошлось (<0.5%)  🟡 близко (<2%)  🔴 не сошлось"))
        self.after(0, self._unlock)

    def _color_converged(self, results):
        """
        Подсвечивает ячейки по степени сходимости.
        Опорное значение — самые мелкие шаги τ и h.
        """
        best_key = (TABLE_TAUS[-1], TABLE_HS[-1])
        best     = results.get(best_key)
        if best is None or abs(best) < 1e-12:
            return

        for (tv, hv), val in results.items():
            rel = abs(val - best) / (abs(best) + 1e-12)
            if rel < 0.005:
                color = GREEN       # сошлось
            elif rel < 0.02:
                color = ACCENT_GOLD # почти сошлось
            else:
                color = ACCENT_HOT  # не сошлось
            self._clbls[(tv, hv)].configure(fg=color)


