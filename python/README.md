# WIP

## Roadmap

### From the v1

* [x] Report IO
* [x] Metagene (Python layer)
* [x] Plots
  * [x] LinePlot (HoloViews / Plotly HTML)
  * [x] HeatMap (HoloViews / Plotly HTML)
  * [x] BoxPlot (HoloViews / Plotly HTML)
* [ ] Gene Body Methylation
* [ ] Converters

### For the v2

* [x] IPC files interface
* [x] IPC indexing interface
* [x] DMR identification algorithm
* [ ] Segmentation algorithm
* [ ] Dimensionality reduction algorithm
* [x] Metagene constructor interface (RegionReader + HcAnnotStore)

TODO: Add get function to AnnotMap

## Metagene Visualization Interface — Design & Usage

Ниже описано, как реализована визуализация метагена поверх backend‑возможностей `bsxplorer2` и как ей пользоваться на реальных `.bsx`/аннотациях.

### Цель

- Создать интерфейс для графического представления результатов анализа `bsxplorer2`.
- Воспользоваться чётким форматом данных `.bsx` и инструментами Rust для эффективной работы с метилированием.
- Восстановить визуализацию метагена из v1: Line plot, Heat map, Box plot, Violin plot, без жёсткой привязки к схеме promoter/body/terminator.

### Подход и Архитектура

- Разделение ответственности и прозрачность компонентов:
  - Данные: `bsx2.plots.data.DiscreteRegionData`, `LinePlotData` — хранение дискретизированных профилей (x/y) и базовые преобразования.
  - Вычисления: `bsx2.plots.metagene.compute_discrete_regions/compute_from_annot` — оркестрация чтения из Rust (`RegionReader`/`HcAnnotStore`), дискретизация посредством `BsxBatch.discretise`.
  - Визуализация: два независимых слоя
    - HoloViews (`line_plot`, `heatmap`, `box_plot`, `violin_plot`) — для интерактивной работы в ноутбуках.
    - Polars + Plotly HTML (`line_html`, `heatmap_html`, `box_html`, `violin_html`) — для генерации standalone HTML.
- Масштабируемость: чтение по регионам стримингом (`RegionReader.iter_contigs`), дискретизация и векторные агрегации (Rust + Polars).
- Гибкость схемы метагена: сегменты произвольны, задаются списком `Segment(name, n_bins)`.

### Работа с Данными Rust (.bsx)

- Источник: `.bsx` — бинарный формат Arrow IPC с индексом, доступ через `RegionReader`.
- Чтение по регионам: `RegionReader.iter_contigs(contigs)` или `RegionReader.query(contig)`. Перед вычислениями можно применить фильтры: `filter_coverage_gt`, `filter_context`, `filter_strand`, `filter_pos_*`.
- Дискретизация внутри региона (Rust): `BsxBatch.discretise(n_bins, AggMethod)` равномерно разбивает геномный отрезок на `n_bins` и агрегирует значения (Mean/Median/Max/Min).
- Получение контуров из аннотаций: `HcAnnotStore.from_gff/.from_bed` + утилита `collect_contigs_from_hcannot`.

### Метаген из Произвольных Участков

- Сегменты задаются пользователем: `segments = [Segment("up", 25), Segment("body", 50), Segment("down", 25)]` или любой иной набор.
- Суммарное число бинов — `sum(s.n_bins)`, `discretise` даёт профили одинаковой длины для всех регионов.
- Ориентация: для ‘-’ цепи профили зеркалируются (опционально), чтобы обеспечить согласованность направления.

### Графики

- HoloViews:
  - `line_plot(reader, contigs=..., segments=..., agg_method=...) -> hv.Curve`
  - `heatmap(...) -> hv.HeatMap`
  - `box_plot(...) -> hv.BoxWhisker`
  - `violin_plot(...) -> hv.Violin`
  - Требуется `hv.extension('matplotlib')` или `'bokeh'` перед использованием `.opts`.
- Plotly HTML (Polars):
  - `line_html(drd, ...) -> str`
  - `heatmap_html(drd, ...) -> str`
  - `box_html(drd, ...) -> str`
  - `violin_html(drd, ...) -> str`

### API Обзор

- Данные: `bsx2.plots.data`
  - `DiscreteRegionData.insert(x: np.ndarray, y: np.ndarray, label: Optional[str])`
  - `DiscreteRegionData.stack_matrix() -> (np.ndarray, list[str])`
  - `LinePlotData.to_curve()` — преобразование в `hv.Curve`.
- Вычисления: `bsx2.plots.metagene`
  - `Segment(name: str, n_bins: int)`, `segments_total_bins`, `segment_boundaries`, `segment_ticks`
  - `compute_discrete_regions(reader, contigs, *, segments, agg_method, ...) -> DiscreteRegionData`
  - `compute_from_annot(reader, annot, *, feature_type, ...) -> DiscreteRegionData`
  - HoloViews: `line_plot`, `heatmap`, `box_plot`, `violin_plot`
- HTML: `bsx2.plots.polars_html`
  - `line_html/heatmap_html/box_html/violin_html`

### Примеры (Кратко)

**HoloViews (RegionReader + Contigs)**

```python
from bsx2.io import RegionReader
from bsx2._bsx2 import Contig, AggMethod
from bsx2.plots import Segment, line_plot
import holoviews as hv; hv.extension('matplotlib')

rr = RegionReader('/path/to/report.bsx')
contigs = [Contig('chr1', 100_000, 120_000, '+')]
segs = [Segment('up', 25), Segment('body', 50), Segment('down', 25)]
curve = line_plot(
    rr,
    contigs=contigs,
    segments=segs,
    agg_method=AggMethod.Mean,
    smooth={"method": "savgol", "window_length": 9, "polyorder": 2},
)
```

**Plotly HTML (из аннотаций)**

```python
from bsx2.io import RegionReader
from bsx2._bsx2 import HcAnnotStore
from bsx2.plots import Segment, compute_from_annot
from bsx2.plots.polars_html import line_html

rr = RegionReader('/path/to/report.bsx')
annot = HcAnnotStore.from_gff('/path/to/annot.gff')
segs = [Segment('up', 25), Segment('body', 50), Segment('down', 25)]
drd = compute_from_annot(rr, annot, segments=segs, limit=100, combine_parts=True)
html = line_html(drd, segments=segs)
open('line.html','w').write(html)
```

### Quick Reference: compute_from_annot

- `reader`: RegionReader над `.bsx`.
- `annot`: HcAnnotStore (from_gff/from_bed).
- `segments`: `list[Segment]`, по умолчанию `[Segment('region', 100)]` (или `up/body/down` при `combine_parts=True`).
- `feature_type`: фильтрация аннотации (например, `'gene'`, применяется при `combine_parts=False`).
- `limit`: ограничение числа регионов (для превью/ускорения).
- `reverse_negative`: инверсия профиля для `'-'` (по умолчанию `True`).
- `labels`: явные подписи регионов (если `None`, используются авто‑лейблы).
- `add_flanks`: добавить фланки gene (up/down) через `annot.add_flanks`.
- `flank_bp`: размер фланков в bp.
- `combine_parts`: собрать BSX1‑стиль (up/body/down) по генам.
- `parts`: порядок частей при `combine_parts=True`.

### Тестирование

- Запуск:

### Используемые Библиотеки и Стиль

- Визуализация: HoloViews (основной), Plotly (HTML без рантайма Python на стороне просмотра).
- Табличные преобразования: Polars (+ PyArrow).
- Валидация типов/структур: возможно подключение `beartype` для рантайм‑валидации (по желанию — `pip install beartype`) и аннотации функций.
- Тесты: `pytest`.
- Управление зависимостями: Poetry (`pyproject.toml`).
- Стиль и форматирование: `ruff` (см. настройки v1), docstrings — стиль NumPy.

### Масштабирование и Рекомендации

- Для больших наборов обязательно использовать `RegionReader.iter_contigs` и фильтры покрытия/контекста для уменьшения NaN и ускорения.
- Для headless‑отчётов предпочтителен Plotly HTML (`*_html*`), для ноутбуков — HoloViews.
- `segments` — свободный выбор структуры, легко адаптировать под любые области интереса (TSS‑центрированные окна, body‑centric и т.д.).

## Python usage (metagene & plots)

Below snippets assume bsx2 is built and importable (see build instructions). The API keeps data prep separate from rendering and supports both HoloViews objects and Plotly HTML.

### Build metagene from RegionReader + Contigs (HoloViews)

```python
from bsx2.io import RegionReader
from bsx2._bsx2 import Contig, AggMethod
from bsx2.plots import Segment, line_plot, heatmap

rr = RegionReader("/path/to/file.bsx")
contigs = [Contig("chr1", 100_000, 120_000, "+"), Contig("chr2", 50_000, 70_000, "-")]
segs = [Segment("up", 25), Segment("body", 50), Segment("down", 25)]

curve = line_plot(rr, contigs=contigs, segments=segs, agg_method=AggMethod.Mean)
hm = heatmap(rr, contigs=contigs, segments=segs, agg_method=AggMethod.Mean)
# curve/hm are HoloViews objects
```

### Build metagene from HcAnnotStore (DiscreteRegionData)

```python
from bsx2.io import RegionReader
from bsx2._bsx2 import HcAnnotStore
from bsx2.plots import Segment, compute_from_annot
from bsx2.plots.polars_html import line_html

rr = RegionReader("/path/to/file.bsx")
annot = HcAnnotStore()  # or load from GFF/BED depending on your pipeline

segs = [Segment("up", 25), Segment("body", 50), Segment("down", 25)]
drd = compute_from_annot(rr, annot, segments=segs, feature_type="gene", combine_parts=True)
html = line_html(drd, segments=segs)
```

### Plotly HTML (standalone)

```python
from bsx2.plots import Segment, compute_from_annot
from bsx2.plots.polars_html import line_html, heatmap_html, box_html, violin_html
from bsx2.io import RegionReader
from bsx2._bsx2 import HcAnnotStore

rr = RegionReader("/path/to/file.bsx")
annot = HcAnnotStore()
segs = [Segment("up", 25), Segment("body", 50), Segment("down", 25)]

drd = compute_from_annot(rr, annot, segments=segs, combine_parts=True)
html_line = line_html(drd, segments=segs)
html_heat = heatmap_html(drd, segments=segs)
html_box = box_html(drd, segments=segs, per_region=True)
html_violin = violin_html(drd, segments=segs, per_region=True)

open("line.html", "w").write(html_line)
open("heatmap.html", "w").write(html_heat)
open("box.html", "w").write(html_box)
open("violin.html", "w").write(html_violin)
```

Notes:
- Segments are arbitrary; you are not constrained to promoter/body/terminator.
- Data prep is independent from visualization; you can take DiscreteRegionData and render via HoloViews or Plotly.
- For large cohorts use Plotly HTML generation (no Python runtime needed to view).

### Quick reference: compute_from_annot parameters

- `reader`: RegionReader over your .bsx file.
- `annot`: HcAnnotStore (e.g. from_gff/from_bed).
- `segments`: list[Segment] — arbitrary segmentation, e.g. [Segment("up",25), Segment("body",50), Segment("down",25)]. If omitted, defaults to [Segment("region", 100)] (or up/body/down when `combine_parts=True`).
- `feature_type`: optional feature filter on annotations (e.g. "gene"; used when `combine_parts=False`).
- `limit`: optional limit on number of regions for faster preview.
- `reverse_negative`: bool — reverse profiles for '-' strand (default True).
- `labels`: optional explicit labels to use for regions.
- `add_flanks`: bool — add upstream/downstream flanks for genes as separate features (default False).
- `flank_bp`: int — flank size in bp for add_flanks (default 2000).
- `combine_parts`: bool — build BSX1-style metagene (up/body/down) per gene.
- `parts`: optional parts order when `combine_parts=True`.

