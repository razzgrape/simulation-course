
import sys
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget,
    QVBoxLayout, QHBoxLayout, QGridLayout,
    QLabel, QSpinBox, QPushButton, QFrame, QSizePolicy,
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal
from PyQt6.QtGui import QFont

from rng_core import get_analysis_results, THEORETICAL_MEAN, THEORETICAL_VARIANCE

SS = """
QMainWindow, QWidget {
    background: #fafafa;
    color: #111;
    font-family: 'Segoe UI';
}
QSpinBox {
    background: #fff;
    border: 1px solid #ddd;
    border-radius: 4px;
    padding: 4px 8px;
    font-size: 13px;
    color: #111;
    min-width: 120px;
}
QSpinBox::up-button, QSpinBox::down-button {
    width: 16px;
    background: #f0f0f0;
    border: none;
}
QPushButton#run {
    background: #111;
    color: #fff;
    border: none;
    border-radius: 4px;
    padding: 8px 24px;
    font-size: 13px;
}
QPushButton#run:hover    { background: #333; }
QPushButton#run:pressed  { background: #000; }
QPushButton#run:disabled { background: #ccc; color: #888; }
"""
class Worker(QThread):
    done = pyqtSignal(dict)

    def __init__(self, n):
        super().__init__()
        self.n = n

    def run(self):
        self.done.emit(get_analysis_results(self.n))


# ── Вспомогательные виджеты ───────────────────────────────────────────────────
def divider():
    f = QFrame()
    f.setFrameShape(QFrame.Shape.HLine)
    f.setStyleSheet("background: #e5e5e5; max-height: 1px; border: none;")
    return f


def label(text, size=12, color="#111", bold=False, mono=False):
    w = QLabel(text)
    f = QFont("Consolas" if mono else "Segoe UI", size)
    f.setBold(bold)
    w.setFont(f)
    w.setStyleSheet(f"color: {color}; background: transparent;")
    return w


class Row(QWidget):
    """Строка: название — значение — отклонение."""
    def __init__(self, name, accent="#111"):
        super().__init__()
        self._accent = accent
        lay = QHBoxLayout(self)
        lay.setContentsMargins(0, 6, 0, 6)
        lay.setSpacing(0)

        self._name  = label(name, size=11, color="#888")
        self._val   = label("—", size=13, color=accent, bold=True, mono=True)
        self._delta = label("—", size=11, color="#aaa", mono=True)

        self._name.setFixedWidth(130)
        self._val.setFixedWidth(120)
        self._delta.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred)

        lay.addWidget(self._name)
        lay.addWidget(self._val)
        lay.addWidget(self._delta)

    def set(self, value: float, delta: float):
        self._val.setText(f"{value:.7f}")
        pct = abs(delta) / THEORETICAL_MEAN * 100
        err_color = "#16a34a" if pct < 0.3 else "#dc2626"
        self._delta.setText(f"Δ {delta:+.2e}  ({pct:.4f}%)")
        self._delta.setStyleSheet(f"color: {err_color}; background: transparent;")


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("RNG")
        self.setFixedSize(560, 480)
        self.setStyleSheet(SS)
        self._build()

    def _build(self):
        root = QWidget()
        self.setCentralWidget(root)
        lay = QVBoxLayout(root)
        lay.setContentsMargins(32, 28, 32, 28)
        lay.setSpacing(0)

        # ── Заголовок
        lay.addWidget(label("Датчик случайных чисел", size=15, bold=True))
        lay.addSpacing(2)
        lay.addWidget(label("Линейный конгруэнтный генератор vs MT19937", size=10, color="#888"))
        lay.addSpacing(20)
        lay.addWidget(divider())
        lay.addSpacing(16)

        # ── Параметры
        p_row = QHBoxLayout()
        p_row.setSpacing(16)
        p_row.addWidget(label("Выборка N", size=11, color="#555"))
        self.spin = QSpinBox()
        self.spin.setRange(1_000, 5_000_000)
        self.spin.setValue(100_000)
        self.spin.setSingleStep(10_000)
        self.spin.setGroupSeparatorShown(True)
        p_row.addWidget(self.spin)
        p_row.addStretch()
        self.btn = QPushButton("Запустить")
        self.btn.setObjectName("run")
        self.btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self.btn.clicked.connect(self._run)
        p_row.addWidget(self.btn)
        lay.addLayout(p_row)
        lay.addSpacing(20)
        lay.addWidget(divider())
        lay.addSpacing(16)

        # ── Шапка таблицы
        hdr = QHBoxLayout()
        hdr.setSpacing(0)
        def hcol(t, w=None):
            l = label(t, size=9, color="#aaa")
            l.setAlignment(Qt.AlignmentFlag.AlignLeft)
            if w: l.setFixedWidth(w)
            return l
        hdr.addWidget(hcol("ПАРАМЕТР", 130))
        hdr.addWidget(hcol("ЗНАЧЕНИЕ", 120))
        hdr.addWidget(hcol("ОТКЛОНЕНИЕ"))
        lay.addLayout(hdr)
        lay.addSpacing(4)
        lay.addWidget(divider())

        # ── Теория
        lay.addSpacing(10)
        lay.addWidget(label("Теоретические", size=9, color="#aaa"))
        lay.addSpacing(4)

        th_grid = QHBoxLayout()
        th_grid.setSpacing(32)
        self._theo_mean = label(f"{THEORETICAL_MEAN:.7f}", size=13, bold=True, mono=True)
        self._theo_var  = label(f"{THEORETICAL_VARIANCE:.7f}", size=13, bold=True, mono=True)
        for ltext, wval in [("Среднее μ", self._theo_mean), ("Дисперсия σ²", self._theo_var)]:
            col = QVBoxLayout()
            col.setSpacing(2)
            col.addWidget(label(ltext, size=9, color="#aaa"))
            col.addWidget(wval)
            th_grid.addLayout(col)
        th_grid.addStretch()
        lay.addLayout(th_grid)
        lay.addSpacing(14)
        lay.addWidget(divider())

        # ── ЛКГ
        lay.addSpacing(10)
        lay.addWidget(label("ЛКГ (собственный)", size=9, color="#2563eb"))
        lay.addSpacing(2)
        self.r_cust_mean = Row("Среднее μ",    "#2563eb")
        self.r_cust_var  = Row("Дисперсия σ²", "#2563eb")
        lay.addWidget(self.r_cust_mean)
        lay.addWidget(self.r_cust_var)
        lay.addSpacing(8)
        lay.addWidget(divider())

        # ── MT19937
        lay.addSpacing(10)
        lay.addWidget(label("MT19937 (random)", size=9, color="#059669"))
        lay.addSpacing(2)
        self.r_built_mean = Row("Среднее μ",    "#059669")
        self.r_built_var  = Row("Дисперсия σ²", "#059669")
        lay.addWidget(self.r_built_mean)
        lay.addWidget(self.r_built_var)

        lay.addStretch()
        lay.addWidget(divider())
        lay.addSpacing(10)

        # ── Статус
        self.status = label("Готов", size=10, color="#aaa")
        lay.addWidget(self.status)

    def _run(self):
        self.btn.setEnabled(False)
        self.status.setText(f"Генерация {self.spin.value():,} значений…")
        self._worker = Worker(self.spin.value())
        self._worker.done.connect(self._done)
        self._worker.start()

    def _done(self, res):
        self.btn.setEnabled(True)
        self.status.setText(f"Готово — N = {res['n']:,}")

        c, b = res["custom"], res["builtin"]
        self.r_cust_mean.set(c["mean"], c["mean"] - THEORETICAL_MEAN)
        self.r_cust_var.set(c["variance"], c["variance"] - THEORETICAL_VARIANCE)
        self.r_built_mean.set(b["mean"], b["mean"] - THEORETICAL_MEAN)
        self.r_built_var.set(b["variance"], b["variance"] - THEORETICAL_VARIANCE)


def main():
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    win = MainWindow()
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()