import sys
import numpy as np
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QLabel, QDoubleSpinBox, QSpinBox, QPushButton, QFormLayout,
    QGroupBox, QSplitter, QGridLayout, QFrame, QTabWidget
)
from PyQt6.QtGui import QFont, QColor
from PyQt6.QtCore import Qt, QTimer
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

# Подключаем твое ядро симулятора
from core import SMOSimulator   

# Меняем шрифт по умолчанию на стандартный системный
APP_FONT = "Segoe UI" if sys.platform == "win32" else "Arial"

class MetricWidget(QFrame):
    """Переименованная карточка статистики в светлом корпоративном стиле"""
    def __init__(self, title, accent_color):
        super().__init__()
        # Вместо темного фиолетового — строгий светло-серый фон с тонкой рамкой
        self.setStyleSheet(f"""
            QFrame {{ 
                background-color: #F8F9FA; 
                border: 1px solid #DEE2E6; 
                border-radius: 6px; 
            }}
        """)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)
        
        self.lbl_title = QLabel(title.upper())
        self.lbl_title.setStyleSheet("color: #6C757D; font-size: 10px; font-weight: bold; letter-spacing: 1px; border: none;")
        self.lbl_title.setAlignment(Qt.AlignmentFlag.AlignLeft)
        
        self.lbl_value = QLabel("0.00")
        self.lbl_value.setStyleSheet(f"color: {accent_color}; font-size: 22px; font-weight: 700; border: none;")
        self.lbl_value.setAlignment(Qt.AlignmentFlag.AlignLeft)
        
        layout.addWidget(self.lbl_title)
        layout.addWidget(self.lbl_value)

    def set_value(self, val):
        self.lbl_value.setText(str(val))

class OperatorStateIndicator(QFrame):
    """Переработанная панель состояния оператора (без смайликов)"""
    def __init__(self):
        super().__init__()
        self.setStyleSheet("background-color: #F8F9FA; border: 1px solid #DEE2E6; border-radius: 6px;")
        layout = QVBoxLayout(self)
        layout.setContentsMargins(10, 10, 10, 10)
        
        self.status_text = QLabel("МОНИТОР СИСТЕМЫ: ОЖИДАНИЕ")
        self.status_text.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.status_text.setStyleSheet("""
            color: #495057; 
            font-size: 13px; 
            font-weight: bold; 
            background-color: #E9ECEF; 
            padding: 8px; 
            border-radius: 4px;
        """)
        layout.addWidget(self.status_text)

    def update_state(self, is_busy):
        # Заменили круглые маркеры на строгую цветовую плашку
        if is_busy:
            self.status_text.setText("СЕРВЕР: ЗАНЯТ")
            self.status_text.setStyleSheet("color: #FFFFFF; background-color: #DC3545; font-size: 13px; font-weight: bold; padding: 8px; border-radius: 4px;")
        else:
            self.status_text.setText("СЕРВЕР: СВОБОДЕН")
            self.status_text.setStyleSheet("color: #FFFFFF; background-color: #28A745; font-size: 13px; font-weight: bold; padding: 8px; border-radius: 4px;")

class SMODashboard(QMainWindow):
    """Главное окно приложения с измененной геометрией панелей"""
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Имитационное моделирование СМО (M/M/1/0)")
        self.resize(1150, 750)
        self.apply_corporate_theme()

        self.sim = None
        self.timer = QTimer()
        self.timer.timeout.connect(self.execution_step)
        self.delay_ms = 40
        self.ticks_per_frame = 4

        self.init_interface()

    def apply_corporate_theme(self):
        # Светлая бизнес-тема вместо неонового "гиковского" Dark Mode
        light_style = f"""
        QMainWindow, QWidget {{ background-color: #FFFFFF; color: #212529; font-family: '{APP_FONT}'; }}
        QGroupBox {{ border: 1px solid #CED4DA; border-radius: 6px; margin-top: 12px; padding-top: 15px; font-weight: bold; background-color: #FFFFFF; }}
        QGroupBox::title {{ subcontrol-origin: margin; subcontrol-position: top left; padding: 0 8px; color: #495057; font-size: 11px; }}
        QSpinBox, QDoubleSpinBox {{ background-color: #FFFFFF; border: 1px solid #CED4DA; padding: 6px; border-radius: 4px; color: #212529; }}
        QSpinBox:focus, QDoubleSpinBox:focus {{ border: 1px solid #0056B3; }}
        
        QPushButton {{ font-weight: bold; font-size: 12px; border-radius: 4px; padding: 10px; color: #FFFFFF; }}
        QPushButton#action_start {{ background-color: #007BFF; }} QPushButton#action_start:hover {{ background-color: #0056B3; }}
        QPushButton#action_pause {{ background-color: #6C757D; }} QPushButton#action_pause:hover {{ background-color: #5A6268; }}
        QPushButton#action_fast {{ background-color: #17A2B8; }} QPushButton#action_fast:hover {{ background-color: #138496; }}
        QPushButton#action_clear {{ background-color: #343A40; }} QPushButton#action_clear:hover {{ background-color: #23272B; }}
        
        QTabWidget::pane {{ border: 1px solid #DEE2E6; border-radius: 6px; background: #FFFFFF; }}
        QTabBar::tab {{ background: #F8F9FA; color: #495057; padding: 8px 16px; border-top-left-radius: 4px; border-top-right-radius: 4px; border: 1px solid #DEE2E6; border-bottom: none; }}
        QTabBar::tab:selected {{ background: #FFFFFF; color: #007BFF; border-top: 2px solid #007BFF; font-weight: bold; }}
        """
        self.setStyleSheet(light_style)

    def init_interface(self):
        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        global_layout = QHBoxLayout(main_widget)
        global_layout.setContentsMargins(12, 12, 12, 12)

        # ---------- ЗЕРКАЛЬНО ПЕРЕВЕРНУЛИ: ЛЕВАЯ ПАНЕЛЬ ТЕПЕРЬ ГРАФИЧЕСКАЯ ----------
        left_side = QWidget()
        left_layout = QVBoxLayout(left_side)
        left_layout.setContentsMargins(0, 0, 0, 0)

        self.operator_monitor = OperatorStateIndicator()
        
        self.view_tabs = QTabWidget()
        self.page_chart = QWidget()
        self.chart_figure = Figure(facecolor='#FFFFFF')
        self.chart_figure.subplots_adjust(bottom=0.15, left=0.1)
        self.chart_canvas = FigureCanvas(self.chart_figure)
        self.chart_axes = self.chart_figure.add_subplot(111)
        
        tmp_layout = QVBoxLayout(self.page_chart)
        tmp_layout.addWidget(self.chart_canvas)
        tmp_layout.setContentsMargins(5, 5, 5, 5)
        self.view_tabs.addTab(self.page_chart, "Вероятности состояний")

        left_layout.addWidget(self.operator_monitor, stretch=1)
        left_layout.addWidget(self.view_tabs, stretch=6)

        # ---------- ПРАВАЯ ПАНЕЛЬ ТЕПЕРЬ УПРАВЛЯЮЩАЯ ----------
        right_side = QWidget()
        right_layout = QVBoxLayout(right_side)
        right_layout.setContentsMargins(0, 0, 0, 0)

        # Блок конфигурации параметров
        config_box = QGroupBox("КОНФИГУРАЦИЯ СМО")
        grid_form = QFormLayout()
        self.input_lambda = QDoubleSpinBox()
        self.input_lambda.setMaximum(9999.0)  # Снимаем ограничение в 99.99
        self.input_lambda.setValue(5.0)
        self.input_lambda.setSingleStep(0.5)

        self.input_mu = QDoubleSpinBox()
        self.input_mu.setMaximum(9999.0)      # Снимаем ограничение в 99.99
        self.input_mu.setValue(6.0)
        self.input_mu.setSingleStep(0.5)
        self.input_total = QSpinBox(); self.input_total.setMaximum(200000); self.input_total.setValue(500)
        grid_form.addRow("Интенсивность потока (λ):", self.input_lambda)
        grid_form.addRow("Интенсивность обработки (μ):", self.input_mu)
        grid_form.addRow("Размер выборки (N):", self.input_total)
        config_box.setLayout(grid_form)

        # Инструменты управления симуляцией
        control_grid = QGridLayout()
        self.btn_run = QPushButton("Запустить"); self.btn_run.setObjectName("action_start")
        self.btn_break = QPushButton("Пауза"); self.btn_break.setObjectName("action_pause")
        self.btn_skip = QPushButton("Экспресс-расчет"); self.btn_skip.setObjectName("action_fast")
        self.btn_wipe = QPushButton("Сбросить"); self.btn_wipe.setObjectName("action_clear")
        
        self.btn_run.clicked.connect(self.start_execution)
        self.btn_break.clicked.connect(self.pause_execution)
        self.btn_skip.clicked.connect(self.run_to_end)
        self.btn_wipe.clicked.connect(self.clear_execution)
        
        control_grid.addWidget(self.btn_run, 0, 0)
        control_grid.addWidget(self.btn_break, 0, 1)
        control_grid.addWidget(self.btn_skip, 1, 0)
        control_grid.addWidget(self.btn_wipe, 1, 1)

        # Мониторинг метрик
        metrics_box = QGroupBox("СВОДНЫЕ ПОКАЗАТЕЛИ")
        metrics_layout = QGridLayout()
        self.m_arrived = MetricWidget("Всего заявок", "#212529")
        self.m_processed = MetricWidget("Успешно обработано", "#28A745")
        self.m_lost = MetricWidget("Отклонено системой", "#DC3545")
        self.m_p_reject = MetricWidget("Коэффициент отказов", "#FD7E14")
        self.m_rho = MetricWidget("Фактическая загрузка", "#6F42C1")
        self.m_avg_time = MetricWidget("Ср. время обработки", "#007BFF")

        metrics_layout.addWidget(self.m_arrived, 0, 0, 1, 2)
        metrics_layout.addWidget(self.m_processed, 1, 0)
        metrics_layout.addWidget(self.m_lost, 1, 1)
        metrics_layout.addWidget(self.m_p_reject, 2, 0)
        metrics_layout.addWidget(self.m_rho, 2, 1)
        metrics_layout.addWidget(self.m_avg_time, 3, 0, 1, 2)
        metrics_box.setLayout(metrics_layout)

        right_layout.addWidget(config_box)
        right_layout.addSpacing(5)
        right_layout.addLayout(control_grid)
        right_layout.addSpacing(5)
        right_layout.addWidget(metrics_box)
        right_layout.addStretch()

        # Объединяем через сплиттер
        window_splitter = QSplitter(Qt.Orientation.Horizontal)
        window_splitter.addWidget(left_side)
        window_splitter.addWidget(right_side)
        window_splitter.setSizes([800, 350])
        global_layout.addWidget(window_splitter)

        self.setup_axis_style(self.chart_axes, "Состояние системы", "Доля времени")

    def setup_axis_style(self, ax, xl, yl):
        ax.set_facecolor('#FFFFFF')
        ax.tick_params(colors='#495057', labelsize=9)
        for border in ['top', 'right']:
            ax.spines[border].set_visible(False)
        for border in ['left', 'bottom']:
            ax.spines[border].set_color('#CED4DA')
        ax.set_xlabel(xl, color='#495057', fontweight='600', family=APP_FONT, fontsize=10)
        ax.set_ylabel(yl, color='#495057', fontweight='600', family=APP_FONT, fontsize=10)
        ax.grid(True, axis='y', linestyle=':', alpha=0.6, color='#6C757D')

    def clear_execution(self):
        self.timer.stop()
        self.sim = SMOSimulator(self.input_lambda.value(), self.input_mu.value(), self.input_total.value())
        self.chart_axes.clear()
        self.setup_axis_style(self.chart_axes, "Состояние системы", "Доля времени")
        self.chart_canvas.draw()
        self.operator_monitor.update_state(False)
        self.refresh_metrics_display()

    def start_execution(self):
        if self.sim is None or self.sim.is_finished:
            self.clear_execution()
        self.timer.start(self.delay_ms)

    def pause_execution(self):
        self.timer.stop()

    def run_to_end(self):
        if self.sim is None:
            self.clear_execution()
        self.timer.stop()
        while not self.sim.is_finished:
            self.sim.step()
        self.render_full_ui(force_graph=True)

    def execution_step(self):
        for _ in range(self.ticks_per_frame):
            if not self.sim.step():
                self.timer.stop()
                self.render_full_ui(force_graph=True)
                return
        self.render_full_ui(force_graph=False)

    def render_full_ui(self, force_graph=True):
        self.operator_monitor.update_state(self.sim.server_busy)
        self.refresh_metrics_display()
        self.generate_plots()

    def generate_plots(self):
        self.chart_axes.clear()
        self.setup_axis_style(self.chart_axes, "Состояние системы", "Доля времени")
        if self.sim and self.sim.current_time > 0:
            total_duration = self.sim.current_time
            keys = sorted(self.sim.state_times.keys())
            frequencies = [self.sim.state_times[k] / total_duration for k in keys]
            
            labels = ['0 (Канал пуст)', '1 (Канал занят)']
            # Строгие цвета графика: синий и серый вместо былой радуги
            palette = ['#6C757D', '#007BFF']
            
            self.chart_axes.bar(keys, frequencies, color=palette, alpha=0.85, width=0.4, align='center')
            self.chart_axes.set_xticks(keys)
            self.chart_axes.set_xticklabels(labels)
            self.chart_axes.set_ylim(0, 1.1)
        self.chart_canvas.draw()

    def refresh_metrics_display(self):
        if self.sim:
            data = self.sim.get_stats()
            self.m_arrived.set_value(str(data['arrived']))
            self.m_processed.set_value(str(data['processed']))
            self.m_lost.set_value(str(data['lost']))

            percentage_lost = data['p_reject'] * 100
            self.m_p_reject.set_value(f"{percentage_lost:.2f} %")

            utilization = data['rho_emp'] * 100
            if data['rho_theor'] >= 1.0:
                self.m_rho.lbl_value.setStyleSheet("color: #DC3545; font-size: 22px; font-weight: 700; border: none;")
            else:
                self.m_rho.lbl_value.setStyleSheet("color: #6F42C1; font-size: 22px; font-weight: 700; border: none;")
            self.m_rho.set_value(f"{utilization:.2f} %")

            self.m_avg_time.set_value(f"{data['avg_service_time']:.4f}")

if __name__ == "__main__":
    app = QApplication(sys.argv)
    app.setStyle('Fusion') # Принудительно ставим единый стиль отображения элементов
    window = SMODashboard()
    window.show()
    sys.exit(app.exec())