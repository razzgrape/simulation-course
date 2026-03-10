import numpy as np
from dataclasses import dataclass, field
from enum import IntEnum
from typing import Optional


class CellState(IntEnum):
    EMPTY = 0
    TREE = 1
    BURNING = 2
    WATER = 3
    BURNT = 4  # промежуточное состояние: выгоревшая клетка


@dataclass
class WindConfig:
    """Конфигурация ветра."""
    enabled: bool = True
    # Направление: 0=N, 1=NE, 2=E, 3=SE, 4=S, 5=SW, 6=W, 7=NW
    direction: int = 2  # по умолчанию восточный
    strength: float = 0.3  # 0..1, дополнительная вероятность воспламенения по ветру


@dataclass
class SimulationConfig:
    """Параметры симуляции."""
    rows: int = 80
    cols: int = 120
    tree_grow_prob: float = 0.01       # p — вероятность роста дерева
    lightning_prob: float = 0.0001     # f — вероятность спонтанного возгорания
    base_spread_prob: float = 0.85     # базовая вероятность распространения огня
    humidity: float = 0.0              # 0..1, снижает вероятность распространения
    wind: WindConfig = field(default_factory=WindConfig)
    # Уклон: матрица высот (None = плоская местность)
    elevation: Optional[np.ndarray] = None
    slope_factor: float = 0.2          # множитель влияния уклона
    burn_duration: int = 2             # сколько шагов клетка горит перед выгоранием

NEIGHBORS = [
    (-1,  0, 0),  # N
    (-1,  1, 1),  # NE
    ( 0,  1, 2),  # E
    ( 1,  1, 3),  # SE
    ( 1,  0, 4),  # S
    ( 1, -1, 5),  # SW
    ( 0, -1, 6),  # W
    (-1, -1, 7),  # NW
]

OPPOSITE_DIR = {0: 4, 1: 5, 2: 6, 3: 7, 4: 0, 5: 1, 6: 2, 7: 3}


class ForestFireSimulation:
    """Двумерный клеточный автомат для моделирования лесных пожаров."""

    def __init__(self, config: SimulationConfig):
        self.config = config
        self.rows = config.rows
        self.cols = config.cols

        self.grid = np.zeros((self.rows, self.cols), dtype=np.int8)
        self.burn_timer = np.zeros((self.rows, self.cols), dtype=np.int8)

        # Генератор высот 
        if config.elevation is None:
            self.elevation = self._generate_elevation()
        else:
            self.elevation = config.elevation

        self.generation = 0
        self.stats_history = []

    def _generate_elevation(self) -> np.ndarray:
        """Генерация карты высот"""
        noise = np.random.randn(self.rows // 8 + 1, self.cols // 8 + 1)
        from scipy.ndimage import zoom
        elevation = zoom(noise, (self.rows / (self.rows // 8 + 1),
                                  self.cols / (self.cols // 8 + 1)))
        elevation = elevation[:self.rows, :self.cols]
        elevation = (elevation - elevation.min()) / (elevation.max() - elevation.min() + 1e-9)
        return elevation

    def initialize_random(self, tree_density: float = 0.6):
        """Случайная инициализация леса"""
        self.grid = np.zeros((self.rows, self.cols), dtype=np.int8)
        self.burn_timer = np.zeros((self.rows, self.cols), dtype=np.int8)
        tree_mask = np.random.random((self.rows, self.cols)) < tree_density
        self.grid[tree_mask] = CellState.TREE
        self.generation = 0
        self.stats_history = []

    def place_water(self, positions: list):
        """Размещение водных преград"""
        for r, c in positions:
            if 0 <= r < self.rows and 0 <= c < self.cols:
                self.grid[r, c] = CellState.WATER

    def add_water_river(self, col: int, width: int = 2):
        """Добавить реку"""
        for r in range(self.rows):
            for dc in range(width):
                c = col + dc + int(2 * np.sin(r / 5.0))  # извилистость
                if 0 <= c < self.cols:
                    self.grid[r, c] = CellState.WATER

    def ignite(self, r: int, c: int):
        """Поджечь клетку."""
        if 0 <= r < self.rows and 0 <= c < self.cols:
            if self.grid[r, c] == CellState.TREE:
                self.grid[r, c] = CellState.BURNING
                self.burn_timer[r, c] = self.config.burn_duration

    def _get_spread_probability(self, from_r: int, from_c: int,
                                  to_r: int, to_c: int, dir_idx: int) -> float:
        """
        Рассчитать вероятность распространения огня.
        Учитывает: базовую вероятность, влажность, ветер, уклон.
        """
        prob = self.config.base_spread_prob

        # Доп правила намбер 1: Влажность 
        # Высокая влажность снижает вероятность
        prob *= (1.0 - self.config.humidity * 0.8)

        # Доп правило намбер 2: Ветер 
        if self.config.wind.enabled and self.config.wind.strength > 0:
            wind_dir = self.config.wind.direction
            # Если направление от горящей к целевой совпадает с ветром — бонус
            if dir_idx == wind_dir:
                prob += self.config.wind.strength
            # Соседние направления тоже получают частичный бонус
            elif dir_idx == (wind_dir + 1) % 8 or dir_idx == (wind_dir - 1) % 8:
                prob += self.config.wind.strength * 0.5
            # Против ветра — штраф
            elif dir_idx == OPPOSITE_DIR[wind_dir]:
                prob -= self.config.wind.strength * 0.6
            elif dir_idx == (OPPOSITE_DIR[wind_dir] + 1) % 8 or dir_idx == (OPPOSITE_DIR[wind_dir] - 1) % 8:
                prob -= self.config.wind.strength * 0.3

        # Доп правило намбер фри: Уклон местности 
        if self.elevation is not None:
            dh = self.elevation[to_r, to_c] - self.elevation[from_r, from_c]
            # Огонь быстрее распространяется вверх по склону
            prob += dh * self.config.slope_factor

        return np.clip(prob, 0.0, 1.0)

    def step(self):
        """Выполнить один шаг автомата."""
        new_grid = self.grid.copy()
        new_burn = self.burn_timer.copy()
        rng = np.random.random((self.rows, self.cols))

        for r in range(self.rows):
            for c in range(self.cols):
                state = self.grid[r, c]

                if state == CellState.BURNING:
                    # Правило намбер уан: горящая клетка уменьшает таймер
                    new_burn[r, c] -= 1
                    if new_burn[r, c] <= 0:
                        new_grid[r, c] = CellState.BURNT

                elif state == CellState.BURNT:
                    # Выгоревшая -> пустая
                    new_grid[r, c] = CellState.EMPTY

                elif state == CellState.TREE:
                    ignited = False
                    # Правило намбер ту: проверяем соседей 
                    for dr, dc, dir_idx in NEIGHBORS:
                        nr, nc = r + dr, c + dc
                        if 0 <= nr < self.rows and 0 <= nc < self.cols:
                            if self.grid[nr, nc] == CellState.BURNING:
                                # Направление распространения: от горящего к дереву
                                spread_dir = OPPOSITE_DIR[dir_idx]
                                prob = self._get_spread_probability(nr, nc, r, c, spread_dir)
                                if rng[r, c] < prob:
                                    new_grid[r, c] = CellState.BURNING
                                    new_burn[r, c] = self.config.burn_duration
                                    ignited = True
                                    break

                    # Правило намбер фри: молния
                    if not ignited:
                        if rng[r, c] < self.config.lightning_prob:
                            new_grid[r, c] = CellState.BURNING
                            new_burn[r, c] = self.config.burn_duration

                elif state == CellState.EMPTY:
                    # Правило намбер фо: рост нового дерева
                    if rng[r, c] < self.config.tree_grow_prob:
                        new_grid[r, c] = CellState.TREE

        self.grid = new_grid
        self.burn_timer = new_burn
        self.generation += 1

        stats = self.get_stats()
        self.stats_history.append(stats)
        return stats

    def get_stats(self) -> dict:
        """Получить текущую статистику."""
        total = self.rows * self.cols
        trees = int(np.sum(self.grid == CellState.TREE))
        burning = int(np.sum(self.grid == CellState.BURNING))
        empty = int(np.sum(self.grid == CellState.EMPTY))
        water = int(np.sum(self.grid == CellState.WATER))
        burnt = int(np.sum(self.grid == CellState.BURNT))
        return {
            'generation': self.generation,
            'trees': trees,
            'burning': burning,
            'empty': empty,
            'water': water,
            'burnt': burnt,
            'total': total,
            'tree_pct': trees / total * 100,
            'burning_pct': burning / total * 100,
        }