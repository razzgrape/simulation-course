import customtkinter as ctk
import tkinter as tk
from tkinter import PhotoImage
import math
import time
import random
import json
import os
from pathlib import Path
from PIL import Image, ImageTk, ImageDraw, ImageOps
from core import get_yes_no, get_magic_8_ball

# ── Тема ──────────────────────────────────────────────────────────────────────
ctk.set_appearance_mode("dark")
ctk.set_default_color_theme("blue")

COLORS = {
    "bg":        "#0D0D1A",
    "surface":   "#13132B",
    "card":      "#1A1A35",
    "accent1":   "#7B5EA7",
    "accent2":   "#E8A838",
    "accent3":   "#4FC3F7",
    "yes":       "#4ADE80",
    "no":        "#F87171",
    "text":      "#F0EAF8",
    "muted":     "#8888AA",
    "border":    "#2A2A50",
}

FONT_MENU  = ("Georgia", 16, "bold")
FONT_SMALL = ("Courier New", 11)

RARITY_COLORS = {
    "common":    "#B0C3D9",
    "uncommon":  "#5E98D9",
    "rare":      "#4B69FF",
    "epic":      "#8847FF",
    "legendary": "#FFD700",
}
RARITY_LABELS = {
    "common":    "Обычный",
    "uncommon":  "Необычный",
    "rare":      "Редкий",
    "epic":      "Эпический",
    "legendary": "Легендарный",
}

ITEM_W  = 130   # ширина слота ленты
ITEM_H  = 130   # высота слота
STRIP_W = 520
WIN_POS = 24    # индекс победителя в ленте


# ══════════════════════════════════════════════════════════════════════════════
#  Загрузка конфига лиц
# ══════════════════════════════════════════════════════════════════════════════

def load_faces(base_dir: str) -> list[dict]:
    """
    Читает faces/config.json — список объектов:
      { "file": "vasya.png", "name": "Вася", "rarity": "common", "weight": 60 }
    Возвращает список с добавленными PIL-изображениями.
    Если config.json нет — сканирует папку и делает всех равновероятными.
    """
    faces_dir = Path(base_dir) / "faces"
    config_path = faces_dir / "config.json"

    if config_path.exists():
        with open(config_path, encoding="utf-8") as f:
            entries = json.load(f)
    else:
        # автоматический режим: все файлы равновероятны
        exts = {".png", ".jpg", ".jpeg", ".webp"}
        files = [p for p in faces_dir.iterdir()
                 if p.suffix.lower() in exts and p.name != "config.json"]
        entries = [
            {"file": p.name, "name": p.stem, "rarity": "common", "weight": 1.0}
            for p in sorted(files)
        ]

    faces = []
    for e in entries:
        img_path = faces_dir / e["file"]
        if not img_path.exists():
            continue
        try:
            img = Image.open(img_path).convert("RGBA")
            img = img.resize((ITEM_W - 8, ITEM_H - 8), Image.LANCZOS)
            faces.append({
                "name":   e.get("name", img_path.stem),
                "rarity": e.get("rarity", "common"),
                "weight": float(e.get("weight", 1.0)),
                "image":  img,
            })
        except Exception:
            pass

    if not faces:
        return []

    total = sum(f["weight"] for f in faces)
    for f in faces:
        f["prob"] = f["weight"] / total

    return faces


def pick_face(faces: list[dict]) -> int:
    """Метод обратного преобразования — выбор лица по весам."""
    alpha = random.random()
    cumulative = 0.0
    for i, f in enumerate(faces):
        cumulative += f["prob"]
        if alpha < cumulative:
            return i
    return len(faces) - 1


def make_placeholder(size: int, color: str = "#2A2A50") -> Image.Image:
    """Заглушка, если лица не загружены."""
    img = Image.new("RGBA", (size, size), color)
    d = ImageDraw.Draw(img)
    d.ellipse([size//4, size//8, 3*size//4, 5*size//8], fill="#4A4A70")
    d.ellipse([size//8, size//2, 7*size//8, size-4], fill="#3A3A60")
    return img


# ══════════════════════════════════════════════════════════════════════════════
#  Вспомогательные виджеты
# ══════════════════════════════════════════════════════════════════════════════

class AnimatedCanvas(tk.Canvas):
    def __init__(self, master, **kw):
        kw.setdefault("bg", COLORS["bg"])
        kw.setdefault("highlightthickness", 0)
        super().__init__(master, **kw)
        self._particles = []
        self._running = False
        self.bind("<Configure>", lambda e: self._init_particles())

    def start(self):
        self._running = True
        self._init_particles()
        self._tick()

    def stop(self):
        self._running = False

    def _init_particles(self):
        self._particles.clear()
        w = self.winfo_width() or 600
        h = self.winfo_height() or 700
        for _ in range(38):
            self._particles.append({
                "x": random.uniform(0, w), "y": random.uniform(0, h),
                "r": random.uniform(1, 3),
                "vx": random.uniform(-0.3, 0.3), "vy": random.uniform(-0.3, 0.3),
                "color": random.choice([COLORS["accent1"], COLORS["accent2"], COLORS["accent3"]]),
            })

    def _tick(self):
        if not self._running:
            return
        self.delete("particle")
        w = self.winfo_width() or 600
        h = self.winfo_height() or 700
        for p in self._particles:
            p["x"] = (p["x"] + p["vx"]) % w
            p["y"] = (p["y"] + p["vy"]) % h
            r = p["r"]
            self.create_oval(p["x"]-r, p["y"]-r, p["x"]+r, p["y"]+r,
                             fill=p["color"], outline="", tags="particle")
        self.after(40, self._tick)


class GlowButton(ctk.CTkButton):
    def __init__(self, master, text, command=None, width=260, height=60,
                 color=COLORS["accent1"], **kw):
        super().__init__(master, text=text, command=command,
                         width=width, height=height, font=FONT_MENU,
                         fg_color=COLORS["card"], hover_color=COLORS["surface"],
                         text_color=COLORS["text"], border_color=color,
                         border_width=2, corner_radius=14)


# ══════════════════════════════════════════════════════════════════════════════
#  Главное меню
# ══════════════════════════════════════════════════════════════════════════════

class MenuScreen(tk.Frame):
    def __init__(self, master, on_yesno, on_magic8, on_faces):
        super().__init__(master, bg=COLORS["bg"])
        self._bg = AnimatedCanvas(self)
        self._bg.place(relwidth=1, relheight=1)

        self._ball_cv = tk.Canvas(self, width=110, height=110,
                                   bg=COLORS["bg"], highlightthickness=0)
        self._ball_cv.place(relx=0.5, rely=0.13, anchor="center")
        self._ang = 0.0
        self._draw_ball()

        tk.Label(self, text="ОРАКУЛ", bg=COLORS["bg"], fg=COLORS["text"],
                 font=("Georgia", 34, "bold")).place(relx=0.5, rely=0.28, anchor="center")
        tk.Label(self, text="Машина вероятностей", bg=COLORS["bg"],
                 fg=COLORS["muted"], font=("Courier New", 12)).place(
                     relx=0.5, rely=0.35, anchor="center")

        GlowButton(self, "⚖  ДА или НЕТ",       command=on_yesno,  color=COLORS["accent3"]).place(relx=0.5, rely=0.49, anchor="center")
        GlowButton(self, "🔮  Шар предсказаний", command=on_magic8, color=COLORS["accent1"]).place(relx=0.5, rely=0.63, anchor="center")
        GlowButton(self, "🎰  Рулетка лиц",      command=on_faces,  color=COLORS["accent2"]).place(relx=0.5, rely=0.77, anchor="center")

        tk.Label(self, text="v1.2  ·  Моделирование случайных событий",
                 bg=COLORS["bg"], fg=COLORS["border"],
                 font=("Courier New", 9)).place(relx=0.5, rely=0.94, anchor="center")

        self._bg.start()
        self._anim_ball()

    def _draw_ball(self):
        c = self._ball_cv
        c.delete("all")
        cx, cy, r = 55, 55, 46
        c.create_oval(cx-r+5, cy-r+8, cx+r+5, cy+r+8, fill="#050510", outline="")
        c.create_oval(cx-r, cy-r, cx+r, cy+r, fill="#1A0A3A",
                      outline=COLORS["accent1"], width=2)
        c.create_oval(cx-18, cy-26, cx+4, cy-8, fill="#FFFFFF",
                      outline="", stipple="gray25")
        off = int(3 * math.sin(self._ang))
        c.create_text(cx+off, cy+off, text="8",
                      fill=COLORS["accent2"], font=("Georgia", 22, "bold"))

    def _anim_ball(self):
        self._ang = (self._ang + 0.05) % (2*math.pi)
        self._draw_ball()
        self.after(40, self._anim_ball)

    def destroy(self):
        self._bg.stop()
        super().destroy()


# ══════════════════════════════════════════════════════════════════════════════
#  Экран «ДА / НЕТ»
# ══════════════════════════════════════════════════════════════════════════════

class YesNoScreen(tk.Frame):
    def __init__(self, master, on_back):
        super().__init__(master, bg=COLORS["bg"])
        self._bg = AnimatedCanvas(self)
        self._bg.place(relwidth=1, relheight=1)

        self._answer_var = tk.StringVar(value="?")
        self._prob_var   = tk.StringVar(value="p = 0.50")
        self._yes = self._no = self._total = 0
        self._animating = False

        tk.Label(self, text="ДА  или  НЕТ", bg=COLORS["bg"], fg=COLORS["text"],
                 font=("Georgia", 28, "bold")).place(relx=0.5, rely=0.07, anchor="center")

        self._ans_lbl = tk.Label(self, textvariable=self._answer_var,
                                  bg=COLORS["bg"], fg=COLORS["accent3"],
                                  font=("Georgia", 72, "bold"))
        self._ans_lbl.place(relx=0.5, rely=0.28, anchor="center")

        tk.Label(self, text="Вероятность «ДА»:", bg=COLORS["bg"],
                 fg=COLORS["muted"], font=FONT_SMALL).place(relx=0.5, rely=0.44, anchor="center")

        self._slider = ctk.CTkSlider(self, from_=0, to=1, number_of_steps=100,
                                      width=280, button_color=COLORS["accent3"],
                                      button_hover_color=COLORS["accent1"],
                                      progress_color=COLORS["accent3"],
                                      fg_color=COLORS["border"],
                                      command=lambda v: self._prob_var.set(f"p = {float(v):.2f}"))
        self._slider.set(0.5)
        self._slider.place(relx=0.5, rely=0.50, anchor="center")

        tk.Label(self, textvariable=self._prob_var, bg=COLORS["bg"],
                 fg=COLORS["accent2"], font=("Courier New", 13, "bold")).place(
                     relx=0.5, rely=0.56, anchor="center")

        GlowButton(self, "  БРОСИТЬ  ", command=self._roll,
                   color=COLORS["yes"], width=220, height=54).place(
                       relx=0.5, rely=0.66, anchor="center")

        sf = tk.Frame(self, bg=COLORS["card"],
                       highlightbackground=COLORS["border"], highlightthickness=1)
        sf.place(relx=0.5, rely=0.82, anchor="center", width=320, height=72)
        self._stat_lbl = tk.Label(sf, text="Статистика появится после первого броска",
                                   bg=COLORS["card"], fg=COLORS["muted"],
                                   font=FONT_SMALL, wraplength=300)
        self._stat_lbl.pack(expand=True)

        GlowButton(self, "← Меню", command=on_back,
                   color=COLORS["muted"], width=130, height=40).place(
                       relx=0.08, rely=0.06, anchor="center")
        self._bg.start()

    def _roll(self):
        if self._animating:
            return
        self._animating = True
        p = self._slider.get()
        chars  = ["?", "!", "ДА", "НЕТ", "?", "ДА", "НЕТ", "?", "!", "?"]
        colors = [COLORS["accent3"], COLORS["accent1"], COLORS["yes"],
                  COLORS["no"], COLORS["accent2"]]

        def step(i):
            if i < len(chars):
                self._answer_var.set(chars[i % len(chars)])
                self._ans_lbl.config(fg=colors[i % len(colors)])
                self.after(60 + i * 15, lambda: step(i+1))
            else:
                result = get_yes_no(p)
                self._answer_var.set(result)
                self._ans_lbl.config(fg=COLORS["yes"] if result == "ДА" else COLORS["no"])
                self._total += 1
                if result == "ДА":
                    self._yes += 1
                else:
                    self._no += 1
                self._stat_lbl.config(
                    text=(f"Бросков: {self._total}   |   "
                          f"ДА: {self._yes} ({self._yes/self._total:.1%})   "
                          f"НЕТ: {self._no} ({self._no/self._total:.1%})"),
                    fg=COLORS["text"])
                self._animating = False
        step(0)

    def destroy(self):
        self._bg.stop()
        super().destroy()


# ══════════════════════════════════════════════════════════════════════════════
#  Экран «Шар предсказаний»
# ══════════════════════════════════════════════════════════════════════════════

class Magic8Screen(tk.Frame):
    def __init__(self, master, on_back):
        super().__init__(master, bg=COLORS["bg"])
        self._bg = AnimatedCanvas(self)
        self._bg.place(relwidth=1, relheight=1)

        self._animating = False
        self._angle     = 0.0
        self._answer    = ""
        self._total     = 0

        tk.Label(self, text="🔮  Шар предсказаний", bg=COLORS["bg"],
                 fg=COLORS["text"], font=("Georgia", 24, "bold")).place(
                     relx=0.5, rely=0.06, anchor="center")

        self._ball_cv = tk.Canvas(self, width=200, height=200,
                                   bg=COLORS["bg"], highlightthickness=0)
        self._ball_cv.place(relx=0.5, rely=0.32, anchor="center")

        tk.Label(self, text="Твой вопрос:", bg=COLORS["bg"],
                 fg=COLORS["muted"], font=FONT_SMALL).place(
                     relx=0.5, rely=0.56, anchor="center")

        self._entry = ctk.CTkEntry(self, width=320, height=38,
                                    placeholder_text="Введи вопрос...",
                                    fg_color=COLORS["card"],
                                    border_color=COLORS["accent1"],
                                    text_color=COLORS["text"],
                                    font=("Courier New", 13))
        self._entry.place(relx=0.5, rely=0.62, anchor="center")
        self._entry.bind("<Return>", lambda e: self._ask())

        GlowButton(self, "✨  СПРОСИТЬ", command=self._ask,
                   color=COLORS["accent1"], width=220, height=52).place(
                       relx=0.5, rely=0.73, anchor="center")

        self._count_lbl = tk.Label(self, text="", bg=COLORS["bg"],
                                    fg=COLORS["muted"], font=FONT_SMALL)
        self._count_lbl.place(relx=0.5, rely=0.84, anchor="center")

        GlowButton(self, "← Меню", command=on_back,
                   color=COLORS["muted"], width=130, height=40).place(
                       relx=0.08, rely=0.06, anchor="center")

        self._bg.start()
        self._draw_ball()
        self._idle_anim()

    def _draw_ball(self, text="8", shake=0):
        c  = self._ball_cv
        cx = 100 + int(shake * math.sin(time.time() * 20))
        cy, r = 100, 88
        c.delete("all")
        c.create_oval(cx-r+10, cy-r+14, cx+r+10, cy+r+14, fill="#050508", outline="")
        for i in range(r, 0, -4):
            d = int(30 + 60 * (1 - i/r))
            col = f"#{d:02x}{max(0,d-15):02x}{min(255,d+40):02x}"
            c.create_oval(cx-i, cy-i, cx+i, cy+i, fill=col, outline="")
        c.create_oval(cx-r, cy-r, cx+r, cy+r, fill="", outline=COLORS["accent1"], width=3)
        ir = 60
        c.create_oval(cx-ir, cy-ir+4, cx+ir, cy+ir+4, fill="#0D1B4A",
                      outline="#1A3080", width=1)
        n = len(text)
        if n <= 1:
            fsize = 30
        elif n <= 6:
            fsize = 22
        elif n <= 14:
            fsize = 14
        elif n <= 22:
            fsize = 12
        else:
            fsize = 10
        c.create_text(cx, cy+4, text=text, fill=COLORS["accent2"],
                      font=("Georgia", fsize, "bold"),
                      width=ir*2-8, justify="center")
        c.create_oval(cx-36, cy-46, cx-10, cy-22, fill="white", outline="", stipple="gray25")
        for deg in range(0, 360, 45):
            a  = math.radians(deg + self._angle * 30)
            px = cx + int((r-5) * math.cos(a))
            py = cy + int((r-5) * math.sin(a))
            c.create_oval(px-2, py-2, px+2, py+2, fill=COLORS["accent1"], outline="")

    def _idle_anim(self):
        if not self._animating:
            self._angle = (self._angle + 0.02) % (2*math.pi)
            self._draw_ball(text=self._answer if self._answer else "8")
        self.after(40, self._idle_anim)

    def _ask(self):
        if self._animating:
            return
        self._animating = True
        self._answer = ""
        self._shake(18, self._reveal)

    def _shake(self, steps, done):
        if steps <= 0:
            self._draw_ball("8", shake=0)
            done()
            return
        self._draw_ball("?", shake=steps * 0.5)
        self.after(40, lambda: self._shake(steps-1, done))

    def _reveal(self):
        answer = get_magic_8_ball()
        self._answer = answer
        self._total += 1
        frames = ["✦", "✧", "✦", answer]

        def show(i):
            if i < len(frames):
                self._draw_ball(frames[i])
                self.after(180, lambda: show(i+1))
            else:
                self._draw_ball(answer)
                self._count_lbl.config(text=f"Вопросов задано: {self._total}")
                self._animating = False
        show(0)

    def destroy(self):
        self._bg.stop()
        super().destroy()


# ══════════════════════════════════════════════════════════════════════════════
#  Экран «Рулетка лиц»
# ══════════════════════════════════════════════════════════════════════════════

class FaceRouletteScreen(tk.Frame):
    def __init__(self, master, on_back, base_dir: str):
        super().__init__(master, bg=COLORS["bg"])

        self._animating  = False
        self._offset     = 0.0
        self._velocity   = 0.0
        self._target_off = 0.0
        self._strip      = []      # список dict из faces
        self._won        = None
        self._history    = []
        self._total      = 0
        self._tk_images  = []      # держим ссылки, чтоб GC не убрал

        # Загружаем лица
        self._faces = load_faces(base_dir)

        # ── Заголовок ──
        tk.Label(self, text="🎰  Рулетка лиц", bg=COLORS["bg"],
                 fg=COLORS["text"], font=("Georgia", 22, "bold")).place(
                     relx=0.5, rely=0.05, anchor="center")

        # ── Лента ──
        wrap = tk.Frame(self, bg=COLORS["accent2"],
                         highlightbackground=COLORS["accent2"],
                         highlightthickness=2)
        wrap.place(relx=0.5, rely=0.22, anchor="center",
                    width=STRIP_W+4, height=ITEM_H+14)

        self._strip_cv = tk.Canvas(wrap, width=STRIP_W, height=ITEM_H+10,
                                    bg=COLORS["surface"], highlightthickness=0)
        self._strip_cv.pack(padx=2, pady=2)

        # ── Маркер ──
        mk = tk.Canvas(self, width=STRIP_W, height=16,
                        bg=COLORS["bg"], highlightthickness=0)
        mk.place(relx=0.5, rely=0.335, anchor="center")
        mid = STRIP_W // 2
        mk.create_polygon(mid-12, 0, mid+12, 0, mid, 14,
                           fill=COLORS["accent2"], outline="")

        # ── Блок результата ──
        rf = tk.Frame(self, bg=COLORS["card"],
                       highlightbackground=COLORS["border"], highlightthickness=1)
        rf.place(relx=0.5, rely=0.45, anchor="center", width=440, height=90)
        self._result_name = tk.Label(rf, text="Нажми «Крутить» — узнай свою судьбу!",
                                      bg=COLORS["card"], fg=COLORS["muted"],
                                      font=("Georgia", 14, "bold"), wraplength=420)
        self._result_name.pack(pady=(10, 0))
        self._result_sub = tk.Label(rf, text="", bg=COLORS["card"],
                                     fg=COLORS["muted"],
                                     font=("Courier New", 10))
        self._result_sub.pack()

        # ── Кнопка ──
        GlowButton(self, "🎲  КРУТИТЬ", command=self._spin_start,
                   color=COLORS["accent2"], width=220, height=52).place(
                       relx=0.5, rely=0.575, anchor="center")

        # ── Вероятности ──
        tk.Label(self, text="Шансы выпадения:", bg=COLORS["bg"],
                 fg=COLORS["muted"], font=("Courier New", 10)).place(
                     relx=0.5, rely=0.655, anchor="center")

        self._prob_cv = tk.Canvas(self, width=460, height=100,
                                   bg=COLORS["bg"], highlightthickness=0)
        self._prob_cv.place(relx=0.5, rely=0.755, anchor="center")
        self._draw_probs()

        # ── История ──
        self._hist_lbl = tk.Label(self, text="", bg=COLORS["bg"],
                                   fg=COLORS["muted"], font=("Courier New", 9),
                                   wraplength=480)
        self._hist_lbl.place(relx=0.5, rely=0.905, anchor="center")

        self._count_lbl = tk.Label(self, text="Крутов: 0", bg=COLORS["bg"],
                                    fg=COLORS["border"], font=("Courier New", 9))
        self._count_lbl.place(relx=0.5, rely=0.955, anchor="center")

        # ── Назад ──
        GlowButton(self, "← Меню", command=on_back,
                   color=COLORS["muted"], width=130, height=40).place(
                       relx=0.08, rely=0.05, anchor="center")

        # Предупреждение, если лица не найдены
        if not self._faces:
            self._result_name.config(
                text="⚠  Папка faces/ не найдена или пуста.\n"
                     "Создай папку faces/ рядом с app.py\n"
                     "и добавь фото + config.json",
                fg=COLORS["no"])

        self._build_idle_strip()
        self._draw_strip(0)

    # ── Таблица вероятностей ──────────────────────────────────────────────────

    def _draw_probs(self):
        c = self._prob_cv
        c.delete("all")
        if not self._faces:
            c.create_text(10, 10, text="Лица не загружены",
                          fill=COLORS["muted"], font=("Courier New", 10), anchor="w")
            return
        x, y = 8, 6
        for f in self._faces:
            col = RARITY_COLORS.get(f["rarity"], COLORS["muted"])
            label = RARITY_LABELS.get(f["rarity"], f["rarity"])
            c.create_text(x, y, text=f"● {f['name']}  —  {f['prob']*100:.1f}%  [{label}]",
                          fill=col, font=("Courier New", 9), anchor="w")
            y += 16
            if y > 95:
                y = 6
                x += 230

    # ── Рендер ленты ─────────────────────────────────────────────────────────

    def _face_to_tk(self, face: dict | None) -> ImageTk.PhotoImage:
        if face is None:
            img = make_placeholder(ITEM_W - 8)
        else:
            img = face["image"].copy()
        tk_img = ImageTk.PhotoImage(img)
        self._tk_images.append(tk_img)   # держим ссылку
        return tk_img

    def _draw_strip(self, offset: float):
        c = self._strip_cv
        c.delete("all")
        self._tk_images.clear()

        center_x = STRIP_W // 2
        first = int(offset // ITEM_W) - 1

        for i in range(first, first + 6):
            if i < 0 or i >= len(self._strip):
                continue
            face = self._strip[i]
            col  = RARITY_COLORS.get(face["rarity"] if face else "common",
                                     COLORS["muted"])
            x    = round(i * ITEM_W - offset + center_x - ITEM_W // 2)

            # Рамка карточки
            c.create_rectangle(x+1, 2, x+ITEM_W-1, ITEM_H+8,
                                fill=COLORS["card"], outline=col, width=2)
            # Полоска редкости снизу
            c.create_rectangle(x+1, ITEM_H+2, x+ITEM_W-1, ITEM_H+8,
                                fill=col, outline="")
            # Фото
            tk_img = self._face_to_tk(face)
            c.create_image(x + ITEM_W//2, (ITEM_H)//2 + 2,
                           image=tk_img, anchor="center")
            # Имя
            c.create_text(x + ITEM_W//2, ITEM_H - 4,
                          text=face["name"] if face else "???",
                          font=("Courier New", 8, "bold"), fill=COLORS["text"],
                          width=ITEM_W - 8, justify="center")

        # Шторки
        c.create_rectangle(0, 0, 60, ITEM_H+10,
                            fill=COLORS["bg"], stipple="gray50", outline="")
        c.create_rectangle(STRIP_W-60, 0, STRIP_W, ITEM_H+10,
                            fill=COLORS["bg"], stipple="gray50", outline="")
        # Центральная линия
        c.create_line(center_x, 0, center_x, ITEM_H+10,
                      fill=COLORS["accent2"], width=1, dash=(4, 4))

    def _build_idle_strip(self):
        if self._faces:
            self._strip = [random.choice(self._faces) for _ in range(8)]
        else:
            self._strip = [{"name": "???", "rarity": "common",
                            "weight": 1, "prob": 1, "image": None}
                           for _ in range(8)]

    def _make_strip(self, winner_idx: int) -> list:
        strip = [random.choice(self._faces) for _ in range(WIN_POS)]
        strip.append(self._faces[winner_idx])
        strip += [random.choice(self._faces) for _ in range(6)]
        return strip

    # ── Кручение ─────────────────────────────────────────────────────────────

    def _spin_start(self):
        if self._animating or not self._faces:
            return
        self._animating = True
        self._result_name.config(text="⏳  Крутим...", fg=COLORS["muted"])
        self._result_sub.config(text="")

        won_idx      = pick_face(self._faces)
        self._won    = self._faces[won_idx]
        self._strip  = self._make_strip(won_idx)

        # Слот i находится под маркером (center_x = STRIP_W//2) когда:
        #   x = i*ITEM_W - offset + center_x - ITEM_W//2 == center_x - ITEM_W//2
        # => offset = i * ITEM_W
        # Добавляем небольшой jitter чтобы победитель останавливался
        # не всегда ровно по центру (реализм CS2-стиля)
        jitter = random.randint(-ITEM_W//3, ITEM_W//3)
        self._target_off = WIN_POS * ITEM_W + jitter
        self._offset     = 0.0
        self._velocity   = 90.0

        self._animate()

    def _animate(self):
        if not self._animating:
            return
        remaining = self._target_off - self._offset
        if remaining <= 0 or self._velocity < 0.35:
            # Снэп: останавливаемся ровно на target, победитель точно под маркером
            self._offset = self._target_off
            self._draw_strip(self._offset)
            self._finish()
            return
        ease = max(0.008, (remaining / max(self._target_off, 1)) ** 0.55)
        self._velocity = max(0.35, 90.0 * ease)
        self._offset  += self._velocity
        self._draw_strip(self._offset)
        self.after(16, self._animate)

    def _finish(self):
        f   = self._won
        col = RARITY_COLORS.get(f["rarity"], COLORS["muted"])
        lbl = RARITY_LABELS.get(f["rarity"], f["rarity"])

        self._result_name.config(text=f"🏆  {f['name']}", fg=col)
        self._result_sub.config(
            text=f"{lbl}  ·  шанс: {f['prob']*100:.2f}%",
            fg=col)

        self._total += 1
        self._count_lbl.config(text=f"Крутов: {self._total}")

        self._history.append(f["name"])
        if len(self._history) > 5:
            self._history.pop(0)
        self._hist_lbl.config(
            text="История: " + "  ·  ".join(self._history))

        self._animating = False


# ══════════════════════════════════════════════════════════════════════════════
#  Главное приложение
# ══════════════════════════════════════════════════════════════════════════════

class App(ctk.CTk):
    def __init__(self):
        super().__init__()
        self.title("Оракул — Машина вероятностей")
        self.geometry("560x700")
        self.minsize(520, 660)
        self.configure(fg_color=COLORS["bg"])
        self._base_dir = os.path.dirname(os.path.abspath(__file__))
        self._current  = None
        self._show_menu()

    def _clear(self):
        if self._current:
            self._current.destroy()
            self._current = None

    def _show_menu(self):
        self._clear()
        f = MenuScreen(self, on_yesno=self._show_yesno,
                       on_magic8=self._show_magic8,
                       on_faces=self._show_faces)
        f.place(relwidth=1, relheight=1)
        self._current = f

    def _show_yesno(self):
        self._clear()
        f = YesNoScreen(self, on_back=self._show_menu)
        f.place(relwidth=1, relheight=1)
        self._current = f

    def _show_magic8(self):
        self._clear()
        f = Magic8Screen(self, on_back=self._show_menu)
        f.place(relwidth=1, relheight=1)
        self._current = f

    def _show_faces(self):
        self._clear()
        f = FaceRouletteScreen(self, on_back=self._show_menu,
                               base_dir=self._base_dir)
        f.place(relwidth=1, relheight=1)
        self._current = f


if __name__ == "__main__":
    App().mainloop()