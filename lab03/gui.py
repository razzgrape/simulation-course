import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
from simulation import (
    ForestFireSimulation, SimulationConfig, WindConfig, CellState
)
from PIL import Image, ImageTk, ImageDraw, ImageFont
import colorsys
import math


# ──────────────────────────────────────────
# Цветовая палитра
# ──────────────────────────────────────────
COLORS = {
    'bg_dark':      '#1a1d23',
    'bg_panel':     '#22262e',
    'bg_card':      '#2a2f38',
    'bg_hover':     '#343a45',
    'accent':       '#ff6b35',
    'accent_light': '#ff8c5a',
    'green':        '#4caf50',
    'green_dark':   '#2e7d32',
    'red':          '#ef5350',
    'yellow':       '#ffc107',
    'blue':         '#42a5f5',
    'blue_dark':    '#1565c0',
    'water':        '#29b6f6',
    'text':         '#e8eaed',
    'text_dim':     '#9aa0a6',
    'text_muted':   '#6b7280',
    'border':       '#3c4049',
    'burnt':        '#5d4037',
}

# Цвета клеток для визуализации
CELL_COLORS = {
    CellState.EMPTY:   (30, 34, 42),       # тёмно-серый
    CellState.TREE:    (46, 125, 50),       # зелёный
    CellState.BURNING: (239, 83, 80),       # красный/оранжевый
    CellState.WATER:   (41, 182, 246),      # голубой
    CellState.BURNT:   (78, 52, 46),        # коричневый
}

# Вариации зелёного для деревьев
TREE_GREENS = [
    (34, 110, 38), (40, 120, 44), (46, 125, 50),
    (56, 135, 60), (38, 115, 42), (50, 130, 54),
    (30, 100, 34), (44, 122, 48), (52, 128, 56),
]

FIRE_COLORS = [
    (255, 87, 34), (239, 83, 80), (255, 152, 0),
    (255, 193, 7), (244, 67, 54), (255, 109, 0),
]

WIND_DIRECTIONS = ['С', 'СВ', 'В', 'ЮВ', 'Ю', 'ЮЗ', 'З', 'СЗ']
WIND_ARROWS = ['↑', '↗', '→', '↘', '↓', '↙', '←', '↖']


class StyledScale(tk.Canvas):
    """Кастомный красивый слайдер."""

    def __init__(self, parent, from_=0, to=1, value=0.5, resolution=0.01,
                 command=None, label="", width=220, **kwargs):
        super().__init__(parent, width=width, height=40,
                        bg=COLORS['bg_card'], highlightthickness=0, **kwargs)
        self.from_ = from_
        self.to = to
        self.value = value
        self.resolution = resolution
        self.command = command
        self.label_text = label
        self.w = width
        self.dragging = False

        self.bind('<Button-1>', self._on_click)
        self.bind('<B1-Motion>', self._on_drag)
        self.bind('<ButtonRelease-1>', self._on_release)
        self._draw()

    def _val_to_x(self, val):
        margin = 12
        return margin + (val - self.from_) / (self.to - self.from_) * (self.w - 2 * margin)

    def _x_to_val(self, x):
        margin = 12
        ratio = (x - margin) / (self.w - 2 * margin)
        ratio = max(0, min(1, ratio))
        val = self.from_ + ratio * (self.to - self.from_)
        if self.resolution >= 1:
            return int(round(val / self.resolution) * self.resolution)
        return round(val / self.resolution) * self.resolution

    def _draw(self):
        self.delete('all')
        margin = 12
        y = 24

        # Трек
        self.create_line(margin, y, self.w - margin, y,
                        fill=COLORS['border'], width=3, capstyle='round')

        # Активная часть
        x = self._val_to_x(self.value)
        self.create_line(margin, y, x, y,
                        fill=COLORS['accent'], width=3, capstyle='round')

        # Кружок
        self.create_oval(x - 7, y - 7, x + 7, y + 7,
                        fill=COLORS['accent'], outline=COLORS['accent_light'], width=2)

        # Подпись и значение
        if isinstance(self.value, float) and self.resolution < 1:
            val_str = f"{self.value:.4f}" if self.resolution < 0.01 else f"{self.value:.2f}"
        else:
            val_str = str(int(self.value))
        self.create_text(margin, 6, text=self.label_text,
                        fill=COLORS['text_dim'], anchor='w', font=('Segoe UI', 8))
        self.create_text(self.w - margin, 6, text=val_str,
                        fill=COLORS['accent_light'], anchor='e', font=('Segoe UI', 8, 'bold'))

    def _on_click(self, event):
        self.dragging = True
        self.value = self._x_to_val(event.x)
        self._draw()
        if self.command:
            self.command(self.value)

    def _on_drag(self, event):
        if self.dragging:
            self.value = self._x_to_val(event.x)
            self._draw()
            if self.command:
                self.command(self.value)

    def _on_release(self, event):
        self.dragging = False

    def set(self, val):
        self.value = val
        self._draw()

    def get(self):
        return self.value


class ForestFireGUI:
    """Основное окно приложения."""

    def __init__(self):
        self.root = tk.Tk()
        self.root.title("🔥 Клеточный автомат — Лесные пожары")
        self.root.configure(bg=COLORS['bg_dark'])
        self.root.geometry("1280x800")
        self.root.minsize(1100, 700)

        # Состояние
        self.running = False
        self.speed = 80  # мс между шагами
        self.cell_size = 6
        self.tool = 'fire'  # fire / tree / water / erase
        self.brush_size = 1

        # Конфигурация и симуляция
        self.config = SimulationConfig(
            rows=100, cols=140,
            tree_grow_prob=0.005,
            lightning_prob=0.00005,
            base_spread_prob=0.7,
            humidity=0.0,
            wind=WindConfig(enabled=True, direction=2, strength=0.2),
            burn_duration=2,
        )
        self.sim = ForestFireSimulation(self.config)
        self.sim.initialize_random(tree_density=0.55)

        self._setup_styles()
        self._build_ui()
        self._render()

    def _setup_styles(self):
        style = ttk.Style()
        style.theme_use('clam')
        style.configure('Dark.TFrame', background=COLORS['bg_dark'])
        style.configure('Panel.TFrame', background=COLORS['bg_panel'])
        style.configure('Card.TFrame', background=COLORS['bg_card'])

    # ──────────────────────────────────────────
    # Построение интерфейса
    # ──────────────────────────────────────────
    def _build_ui(self):
        # Главный контейнер
        main = tk.Frame(self.root, bg=COLORS['bg_dark'])
        main.pack(fill='both', expand=True, padx=8, pady=8)

        # Левая панель управления
        self._build_control_panel(main)

        # Центральная область
        center = tk.Frame(main, bg=COLORS['bg_dark'])
        center.pack(side='left', fill='both', expand=True, padx=(8, 0))

        # Верхняя панель инструментов
        self._build_toolbar(center)

        # Канвас
        canvas_frame = tk.Frame(center, bg=COLORS['border'], padx=2, pady=2)
        canvas_frame.pack(fill='both', expand=True, pady=(6, 0))

        self.canvas = tk.Canvas(canvas_frame, bg=COLORS['bg_dark'],
                                highlightthickness=0, cursor='crosshair')
        self.canvas.pack(fill='both', expand=True)
        self.canvas.bind('<Button-1>', self._on_canvas_click)
        self.canvas.bind('<B1-Motion>', self._on_canvas_drag)
        self.canvas.bind('<Configure>', self._on_canvas_resize)

        # Статусбар
        self._build_statusbar(center)

    def _build_control_panel(self, parent):
        panel = tk.Frame(parent, bg=COLORS['bg_panel'], width=260)
        panel.pack(side='left', fill='y')
        panel.pack_propagate(False)

        # Заголовок
        title_frame = tk.Frame(panel, bg=COLORS['bg_panel'])
        title_frame.pack(fill='x', padx=12, pady=(14, 4))
        tk.Label(title_frame, text="🌲 ЛЕСНЫЕ ПОЖАРЫ",
                font=('Segoe UI', 13, 'bold'), fg=COLORS['accent'],
                bg=COLORS['bg_panel']).pack(anchor='w')
        tk.Label(title_frame, text="Клеточный автомат",
                font=('Segoe UI', 9), fg=COLORS['text_dim'],
                bg=COLORS['bg_panel']).pack(anchor='w')

        # Разделитель
        self._separator(panel)

        # ── Управление симуляцией ──
        self._section_label(panel, "⏱  УПРАВЛЕНИЕ")

        btn_frame = tk.Frame(panel, bg=COLORS['bg_panel'])
        btn_frame.pack(fill='x', padx=12, pady=4)

        self.btn_play = self._button(btn_frame, "▶ Старт", self._toggle_play,
                                     bg=COLORS['green_dark'])
        self.btn_play.pack(side='left', expand=True, fill='x', padx=(0, 3))

        self._button(btn_frame, "⏭ Шаг", self._step_once).pack(
            side='left', expand=True, fill='x', padx=(3, 3))

        self._button(btn_frame, "🔄", self._reset, width=3).pack(
            side='left', padx=(3, 0))

        # Скорость
        self.speed_slider = StyledScale(
            panel, from_=10, to=300, value=self.speed,
            resolution=10, command=self._on_speed_change, label="Скорость (мс)")
        self.speed_slider.pack(padx=12, pady=(6, 2))

        self._separator(panel)

        # ── Параметры ──
        self._section_label(panel, "⚙  ПАРАМЕТРЫ")

        # Scrollable parameters
        params_canvas = tk.Canvas(panel, bg=COLORS['bg_panel'],
                                  highlightthickness=0, width=236)
        scrollbar = tk.Scrollbar(panel, orient='vertical', command=params_canvas.yview)
        params_inner = tk.Frame(params_canvas, bg=COLORS['bg_panel'])

        params_inner.bind('<Configure>',
                         lambda e: params_canvas.configure(scrollregion=params_canvas.bbox('all')))
        params_canvas.create_window((0, 0), window=params_inner, anchor='nw')
        params_canvas.configure(yscrollcommand=scrollbar.set)

        params_canvas.pack(side='left', fill='both', expand=True)
        scrollbar.pack(side='right', fill='y')

        # Bind mousewheel to params canvas
        def _on_mousewheel(event):
            params_canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        params_canvas.bind_all('<MouseWheel>', _on_mousewheel)

        self._build_params(params_inner)

    def _build_params(self, parent):
        # Рост деревьев
        self.grow_slider = StyledScale(
            parent, from_=0, to=0.05, value=self.config.tree_grow_prob,
            resolution=0.001, command=lambda v: setattr(self.config, 'tree_grow_prob', v),
            label="Рост деревьев (p)")
        self.grow_slider.pack(padx=8, pady=4)

        # Молния
        self.lightning_slider = StyledScale(
            parent, from_=0, to=0.001, value=self.config.lightning_prob,
            resolution=0.00001, command=lambda v: setattr(self.config, 'lightning_prob', v),
            label="Молния (f)")
        self.lightning_slider.pack(padx=8, pady=4)

        # Вероятность распространения
        self.spread_slider = StyledScale(
            parent, from_=0, to=1, value=self.config.base_spread_prob,
            resolution=0.05, command=lambda v: setattr(self.config, 'base_spread_prob', v),
            label="Распространение огня")
        self.spread_slider.pack(padx=8, pady=4)

        # Длительность горения
        self.burn_slider = StyledScale(
            parent, from_=1, to=5, value=self.config.burn_duration,
            resolution=1, command=lambda v: setattr(self.config, 'burn_duration', int(v)),
            label="Длительность горения")
        self.burn_slider.pack(padx=8, pady=4)

        # Плотность леса (при сбросе)
        self.density_slider = StyledScale(
            parent, from_=0, to=1, value=0.55,
            resolution=0.05, label="Плотность леса")
        self.density_slider.pack(padx=8, pady=4)

        self._mini_separator(parent)

        # ── Ветер ──
        self._section_label(parent, "💨  ВЕТЕР", padx=8)

        wind_frame = tk.Frame(parent, bg=COLORS['bg_panel'])
        wind_frame.pack(fill='x', padx=8, pady=4)

        self.wind_enabled_var = tk.BooleanVar(value=self.config.wind.enabled)
        cb = tk.Checkbutton(wind_frame, text="Включён", variable=self.wind_enabled_var,
                           command=self._on_wind_toggle,
                           fg=COLORS['text'], bg=COLORS['bg_panel'],
                           selectcolor=COLORS['bg_card'], activebackground=COLORS['bg_panel'],
                           font=('Segoe UI', 9))
        cb.pack(anchor='w')

        # Компас для направления ветра
        self._build_wind_compass(parent)

        self.wind_strength_slider = StyledScale(
            parent, from_=0, to=0.8, value=self.config.wind.strength,
            resolution=0.05, command=lambda v: setattr(self.config.wind, 'strength', v),
            label="Сила ветра")
        self.wind_strength_slider.pack(padx=8, pady=4)

        self._mini_separator(parent)

        # ── Влажность ──
        self._section_label(parent, "💧  ВЛАЖНОСТЬ", padx=8)

        self.humidity_slider = StyledScale(
            parent, from_=0, to=1, value=self.config.humidity,
            resolution=0.05, command=lambda v: setattr(self.config, 'humidity', v),
            label="Уровень влажности")
        self.humidity_slider.pack(padx=8, pady=4)

        self._mini_separator(parent)

        # ── Уклон ──
        self._section_label(parent, "⛰  УКЛОН МЕСТНОСТИ", padx=8)

        self.slope_slider = StyledScale(
            parent, from_=0, to=1, value=self.config.slope_factor,
            resolution=0.05, command=lambda v: setattr(self.config, 'slope_factor', v),
            label="Влияние уклона")
        self.slope_slider.pack(padx=8, pady=4)

        regen_btn = self._button(parent, "🔄 Новый рельеф", self._regen_elevation)
        regen_btn.pack(padx=8, pady=(2, 8), fill='x')

    def _build_wind_compass(self, parent):
        """Мини-компас для выбора направления ветра."""
        compass = tk.Canvas(parent, width=120, height=120,
                           bg=COLORS['bg_card'], highlightthickness=0)
        compass.pack(padx=8, pady=4)
        self.compass_canvas = compass

        self._draw_compass()

    def _draw_compass(self):
        c = self.compass_canvas
        c.delete('all')
        cx, cy, r = 60, 60, 42

        # Круг фона
        c.create_oval(cx-r, cy-r, cx+r, cy+r,
                      fill=COLORS['bg_dark'], outline=COLORS['border'], width=2)

        # Направления
        for i, (label, arrow) in enumerate(zip(WIND_DIRECTIONS, WIND_ARROWS)):
            angle = -90 + i * 45  # начинаем с севера
            rad = math.radians(angle)
            x = cx + (r - 8) * math.cos(rad)
            y = cy + (r - 8) * math.sin(rad)

            is_selected = (i == self.config.wind.direction)
            color = COLORS['accent'] if is_selected else COLORS['text_muted']
            font_weight = 'bold' if is_selected else 'normal'

            tag = f'dir_{i}'
            c.create_text(x, y, text=arrow, fill=color,
                         font=('Segoe UI', 11, font_weight), tags=tag)
            c.tag_bind(tag, '<Button-1>',
                      lambda e, idx=i: self._set_wind_direction(idx))

        # Центр — текущее направление
        wd = self.config.wind.direction
        c.create_text(cx, cy, text=WIND_ARROWS[wd],
                     fill=COLORS['accent_light'],
                     font=('Segoe UI', 16, 'bold'))

    def _set_wind_direction(self, idx):
        self.config.wind.direction = idx
        self._draw_compass()

    def _build_toolbar(self, parent):
        bar = tk.Frame(parent, bg=COLORS['bg_panel'], height=50)
        bar.pack(fill='x')
        bar.pack_propagate(False)

        # Инструменты рисования
        tk.Label(bar, text="Инструменты:", font=('Helvetica', 10, 'bold'),
                fg=COLORS['text_dim'], bg=COLORS['bg_panel']).pack(
                    side='left', padx=(10, 6))

        tools = [
            ('🔥 Огонь', 'fire'),
            ('🌲 Дерево', 'tree'),
            ('💧 Вода', 'water'),
            ('⬜ Стереть', 'erase'),
        ]
        self.tool_buttons = {}
        for text, tool_id in tools:
            btn = self._tool_button(bar, text, tool_id)
            btn.pack(side='left', padx=2, pady=6)
            self.tool_buttons[tool_id] = btn

        self._highlight_tool('fire')

        # Размер кисти
        tk.Label(bar, text="  Кисть:", font=('Helvetica', 10, 'bold'),
                fg=COLORS['text_dim'], bg=COLORS['bg_panel']).pack(side='left', padx=(10, 4))
        self.brush_buttons = {}
        for size in [1, 3, 5]:
            bg_color = '#3a4050'
            bf = tk.Frame(bar, bg=bg_color, padx=1, pady=1,
                         highlightbackground='#555', highlightthickness=1)
            bl = tk.Label(bf, text=str(size), width=3,
                         font=('Helvetica', 11, 'bold'),
                         fg='#ffffff', bg=bg_color,
                         padx=4, pady=2, cursor='hand2')
            bl.pack()
            bl.bind('<Button-1>', lambda e, s=size: self._set_brush(s))
            bf.pack(side='left', padx=2, pady=6)
            bf._label = bl
            self.brush_buttons[size] = bf
        self._highlight_brush(1)

        # Добавить реку
        self._button(bar, "〰 Река", self._add_river, width=8).pack(
            side='left', padx=(15, 2), pady=6)

        # Статистика справа
        self.stats_label = tk.Label(bar, text="",
                                    font=('Consolas', 9),
                                    fg=COLORS['text_dim'],
                                    bg=COLORS['bg_panel'])
        self.stats_label.pack(side='right', padx=10)

    def _build_statusbar(self, parent):
        bar = tk.Frame(parent, bg=COLORS['bg_panel'], height=32)
        bar.pack(fill='x', pady=(4, 0))
        bar.pack_propagate(False)

        self.gen_label = tk.Label(bar, text="Поколение: 0",
                                 font=('Consolas', 9, 'bold'),
                                 fg=COLORS['accent_light'],
                                 bg=COLORS['bg_panel'])
        self.gen_label.pack(side='left', padx=10)

        # Легенда
        legend_items = [
            ("Пусто", COLORS['bg_dark']),
            ("Дерево", '#2e7d32'),
            ("Горит", '#ef5350'),
            ("Вода", '#29b6f6'),
            ("Выгорело", '#5d4037'),
        ]
        for text, color in legend_items:
            f = tk.Frame(bar, bg=COLORS['bg_panel'])
            f.pack(side='left', padx=(8, 0))
            tk.Canvas(f, width=10, height=10, bg=color,
                     highlightthickness=1, highlightbackground=COLORS['border']
                     ).pack(side='left', padx=(0, 3))
            tk.Label(f, text=text, font=('Segoe UI', 8),
                    fg=COLORS['text_dim'], bg=COLORS['bg_panel']).pack(side='left')

        # Доп. правила справа
        rules_text = "Доп. правила: Ветер | Влажность | Преграды | Уклон"
        tk.Label(bar, text=rules_text, font=('Segoe UI', 8),
                fg=COLORS['text_muted'], bg=COLORS['bg_panel']).pack(
                    side='right', padx=10)

    # ──────────────────────────────────────────
    # UI-помощники
    # ──────────────────────────────────────────
    def _button(self, parent, text, command, bg=None, width=None, **kwargs):
        """Label-based button — macOS tk.Button ignores bg/fg colors."""
        bg_color = bg or '#3a4050'
        frame = tk.Frame(parent, bg=bg_color, padx=1, pady=1,
                        highlightbackground='#555', highlightthickness=1)
        lbl = tk.Label(frame, text=text,
                      font=('Helvetica', 12, 'bold'),
                      fg='#ffffff', bg=bg_color,
                      padx=12, pady=5, cursor='hand2')
        lbl.pack(fill='both', expand=True)
        if width:
            lbl.configure(width=width)
        # Hover эффект
        hover_color = '#' + ''.join(
            f'{min(255, int(bg_color[i:i+2], 16) + 25):02x}'
            for i in (1, 3, 5)
        )
        lbl.bind('<Enter>', lambda e: lbl.configure(bg=hover_color))
        lbl.bind('<Leave>', lambda e: lbl.configure(bg=bg_color))
        lbl.bind('<Button-1>', lambda e: command())
        frame._label = lbl
        frame._bg_color = bg_color
        # Метод config для совместимости с btn_play.config(...)
        _orig_config = frame.configure
        def _config_wrapper(**kw):
            if 'text' in kw:
                lbl.configure(text=kw.pop('text'))
            if 'bg' in kw:
                new_bg = kw.pop('bg')
                frame.configure(bg=new_bg)
                lbl.configure(bg=new_bg)
                frame._bg_color = new_bg
                nonlocal hover_color
                hover_color = '#' + ''.join(
                    f'{min(255, int(new_bg[i:i+2], 16) + 25):02x}'
                    for i in (1, 3, 5)
                )
                lbl.bind('<Enter>', lambda e: lbl.configure(bg=hover_color))
                lbl.bind('<Leave>', lambda e: lbl.configure(bg=new_bg))
            if 'fg' in kw:
                lbl.configure(fg=kw.pop('fg'))
            if kw:
                _orig_config(**kw)
        frame.config = _config_wrapper
        return frame

    def _tool_button(self, parent, text, tool_id):
        bg_color = '#3a4050'
        frame = tk.Frame(parent, bg=bg_color, padx=1, pady=1,
                        highlightbackground='#555', highlightthickness=1)
        lbl = tk.Label(frame, text=text,
                      font=('Helvetica', 11, 'bold'),
                      fg='#ffffff', bg=bg_color,
                      padx=10, pady=3, cursor='hand2')
        lbl.pack(fill='both', expand=True)
        lbl.bind('<Button-1>', lambda e: self._set_tool(tool_id))
        frame._label = lbl
        frame._bg_color = bg_color
        return frame

    def _highlight_tool(self, tool_id):
        for tid, frame in self.tool_buttons.items():
            lbl = frame._label
            if tid == tool_id:
                frame.configure(bg=COLORS['accent'],
                               highlightbackground=COLORS['accent_light'])
                lbl.configure(bg=COLORS['accent'], fg='#ffffff')
            else:
                frame.configure(bg='#3a4050', highlightbackground='#555')
                lbl.configure(bg='#3a4050', fg='#cccccc')

    def _section_label(self, parent, text, padx=12):
        tk.Label(parent, text=text, font=('Segoe UI', 9, 'bold'),
                fg=COLORS['text_dim'], bg=COLORS['bg_panel']).pack(
                    anchor='w', padx=padx, pady=(8, 2))

    def _separator(self, parent):
        tk.Frame(parent, bg=COLORS['border'], height=1).pack(
            fill='x', padx=12, pady=8)

    def _mini_separator(self, parent):
        tk.Frame(parent, bg=COLORS['border'], height=1).pack(
            fill='x', padx=20, pady=6)

    # ──────────────────────────────────────────
    # Рендеринг
    # ──────────────────────────────────────────
    def _render(self):
        """Отрисовка сетки как изображения на canvas."""
        grid = self.sim.grid
        rows, cols = grid.shape

        # Создаём numpy-массив цветов
        img_array = np.zeros((rows, cols, 3), dtype=np.uint8)

        for state, color in CELL_COLORS.items():
            mask = (grid == state)
            img_array[mask] = color

        # Вариации для деревьев (используем хеш позиции)
        tree_mask = (grid == CellState.TREE)
        tree_indices = np.argwhere(tree_mask)
        for r, c in tree_indices:
            idx = (r * 7 + c * 13) % len(TREE_GREENS)
            img_array[r, c] = TREE_GREENS[idx]

        # Мерцание огня
        fire_mask = (grid == CellState.BURNING)
        fire_indices = np.argwhere(fire_mask)
        for r, c in fire_indices:
            idx = (r + c + self.sim.generation) % len(FIRE_COLORS)
            img_array[r, c] = FIRE_COLORS[idx]

        # Наложение карты высот (подсветка) — если уклон включен
        if self.config.slope_factor > 0 and self.sim.elevation is not None:
            elev = self.sim.elevation
            brightness = 0.7 + 0.3 * elev  # от 0.7 до 1.0
            for ch in range(3):
                img_array[:, :, ch] = np.clip(
                    img_array[:, :, ch] * brightness, 0, 255
                ).astype(np.uint8)

        # Масштабирование
        canvas_w = self.canvas.winfo_width()
        canvas_h = self.canvas.winfo_height()
        if canvas_w < 10 or canvas_h < 10:
            return

        img = Image.fromarray(img_array, 'RGB')
        # Вычисляем размер ячейки чтобы вписать в канвас
        cell_w = max(1, canvas_w // cols)
        cell_h = max(1, canvas_h // rows)
        self.cell_size = min(cell_w, cell_h)

        new_w = cols * self.cell_size
        new_h = rows * self.cell_size
        img = img.resize((new_w, new_h), Image.NEAREST)

        self._photo = ImageTk.PhotoImage(img)
        self.canvas.delete('all')
        self.canvas.create_image(canvas_w // 2, canvas_h // 2,
                                image=self._photo, anchor='center')

        # Обновить статистику
        self._update_stats()

    def _update_stats(self):
        stats = self.sim.get_stats()
        self.gen_label.config(text=f"Поколение: {stats['generation']}")
        self.stats_label.config(
            text=f"🌲 {stats['trees']}  ({stats['tree_pct']:.1f}%)  "
                 f"🔥 {stats['burning']}  ({stats['burning_pct']:.1f}%)  "
                 f"💧 {stats['water']}"
        )

    # ──────────────────────────────────────────
    # Обработчики событий
    # ──────────────────────────────────────────
    def _toggle_play(self):
        self.running = not self.running
        if self.running:
            self.btn_play.config(text="⏸ Пауза", bg='#d32f2f',
                                fg='#ffffff')
            self._run_loop()
        else:
            self.btn_play.config(text="▶ Старт", bg='#2e7d32',
                                fg='#ffffff')

    def _run_loop(self):
        if self.running:
            self.sim.step()
            self._render()
            self.root.after(self.speed, self._run_loop)

    def _step_once(self):
        self.sim.step()
        self._render()

    def _reset(self):
        self.running = False
        self.btn_play.config(text="▶ Старт", bg='#2e7d32', fg='#ffffff')
        density = self.density_slider.get()
        self.sim = ForestFireSimulation(self.config)
        self.sim.initialize_random(tree_density=density)
        self._render()

    def _on_speed_change(self, val):
        self.speed = int(val)

    def _set_tool(self, tool_id):
        self.tool = tool_id
        self._highlight_tool(tool_id)

    def _set_brush(self, size):
        self.brush_size = size
        self._highlight_brush(size)

    def _highlight_brush(self, size):
        for s, frame in self.brush_buttons.items():
            lbl = frame._label
            if s == size:
                frame.configure(bg=COLORS['accent'],
                               highlightbackground=COLORS['accent_light'])
                lbl.configure(bg=COLORS['accent'], fg='#ffffff')
            else:
                frame.configure(bg='#3a4050', highlightbackground='#555')
                lbl.configure(bg='#3a4050', fg='#cccccc')

    def _canvas_to_grid(self, event):
        """Преобразование координат канваса в координаты сетки."""
        canvas_w = self.canvas.winfo_width()
        canvas_h = self.canvas.winfo_height()
        rows, cols = self.sim.rows, self.sim.cols
        grid_w = cols * self.cell_size
        grid_h = rows * self.cell_size
        offset_x = (canvas_w - grid_w) / 2
        offset_y = (canvas_h - grid_h) / 2

        c = int((event.x - offset_x) / self.cell_size)
        r = int((event.y - offset_y) / self.cell_size)
        return r, c

    def _apply_tool(self, r, c):
        half = self.brush_size // 2
        for dr in range(-half, half + 1):
            for dc in range(-half, half + 1):
                rr, cc = r + dr, c + dc
                if 0 <= rr < self.sim.rows and 0 <= cc < self.sim.cols:
                    if self.tool == 'fire':
                        self.sim.ignite(rr, cc)
                    elif self.tool == 'tree':
                        if self.sim.grid[rr, cc] != CellState.WATER:
                            self.sim.grid[rr, cc] = CellState.TREE
                    elif self.tool == 'water':
                        self.sim.grid[rr, cc] = CellState.WATER
                    elif self.tool == 'erase':
                        self.sim.grid[rr, cc] = CellState.EMPTY

    def _on_canvas_click(self, event):
        r, c = self._canvas_to_grid(event)
        self._apply_tool(r, c)
        self._render()

    def _on_canvas_drag(self, event):
        r, c = self._canvas_to_grid(event)
        self._apply_tool(r, c)
        self._render()

    def _on_canvas_resize(self, event):
        self._render()

    def _on_wind_toggle(self):
        self.config.wind.enabled = self.wind_enabled_var.get()

    def _add_river(self):
        col = np.random.randint(20, self.sim.cols - 20)
        self.sim.add_water_river(col, width=3)
        self._render()

    def _regen_elevation(self):
        self.sim.elevation = self.sim._generate_elevation()
        self._render()

    # ──────────────────────────────────────────
    # Запуск
    # ──────────────────────────────────────────
    def run(self):
        self.root.mainloop()


if __name__ == '__main__':
    app = ForestFireGUI()
    app.run()