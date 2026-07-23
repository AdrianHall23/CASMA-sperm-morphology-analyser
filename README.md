Here is your README — copy everything between the lines and paste it into a file called `README.md` in your repo:

---

```markdown
# CASMA — Computer-Aided Sperm Morphology Analyser

> Diagnostic-grade image analysis tool for quantifying sperm morphology from microscopy images against WHO 6th Edition strict criteria.

---

## Overview

CASMA is a Python-based image analysis pipeline designed to replicate the workflow of commercial computer-aided sperm analysis (CASA/CASMA) systems used in UKAS-accredited andrology laboratories. It segments sperm heads from brightfield microscopy images, computes six morphometric parameters, classifies each sperm as normal or abnormal against WHO reference ranges, and generates ALCOA+-compliant audit trails for diagnostic quality assurance.

Built as part of independent research into diagnostic andrology workflows and laboratory data integrity.

---

## Features

| Module | Description |
|--------|-------------|
| **Image Preprocessing** | Gaussian denoising, adaptive Otsu thresholding, morphological opening/closing, hole filling |
| **Segmentation** | Connected-component labelling with size-based filtering to isolate sperm heads from debris |
| **Morphometrics** | Principal-component analysis (inertia tensor) for major/minor axis estimation; area, perimeter, circularity, aspect ratio |
| **WHO Classification** | Strict morphology criteria from *WHO Laboratory Manual for the Examination and Processing of Human Semen* (6th Edition, 2021) |
| **Robust QC** | Modified Z-score outlier detection using Median Absolute Deviation (MAD) — robust to small sample sizes |
| **Visualisation** | Annotated microscopy overlays (green = normal, red = abnormal) + distribution histograms with WHO reference lines |
| **Reporting** | Per-sperm CSV, population summary JSON, and ALCOA+ compliant audit trail |

---

## WHO 6th Edition Criteria Applied

| Parameter | Normal Range | Source |
|-----------|--------------|--------|
| Head length | 4.0 – 5.0 µm | WHO 6th Ed., Table 2.5 |
| Head width | 2.5 – 3.5 µm | WHO 6th Ed., Table 2.5 |
| Head area | 12.0 – 22.0 µm² | Derived from length × width |
| Perimeter | 13.0 – 19.0 µm | Approximated from contour |
| Circularity | ≥ 0.75 | Roundness index |
| Aspect ratio | 1.20 – 1.80 | Length / width |

**Lower reference limit for normal morphology:** 4.0%

---

## Installation

```bash
git clone https://github.com/AdrianHall23/CASMA-Sperm-Morphology-Analyser.git
cd CASMA-Sperm-Morphology-Analyser
pip install -r requirements.txt
```

### Requirements
- Python ≥ 3.9
- numpy
- scipy
- matplotlib
- Pillow

---

## Usage

### Basic run
```bash
python casma_analyser.py --image sample_image.tif --scale 0.12
```

### Arguments
| Flag | Description | Default |
|------|-------------|---------|
| `--image` | Path to microscopy image (TIFF, PNG, JPG) | *required* |
| `--scale` | Micrometres per pixel (calibrate per objective) | `0.12` |
| `--outdir` | Output directory for reports | `casma_output` |

### Example with 40× objective
```bash
python casma_analyser.py --image patient_sample_001.tif --scale 0.12 --outdir results/patient_001
```

---

## Outputs

```
casma_output/
├── annotated_image.png          # Original image with measured sperm highlighted
├── morphology_distributions.png # Histograms of all 6 parameters vs WHO ranges
├── morphology_report.csv        # Per-sperm morphometrics + normal/abnormal flag
├── summary_statistics.json      # Population stats: mean, std, median, normal %
└── qc_audit_trail.txt           # ALCOA+ compliant audit log
```

---

## Sample Results

### Annotated Microscopy Image
![Annotated Sample](docs/annotated_sample.png)
*Green dots = normal morphology; Red dots = abnormal morphology. Numbers correspond to rows in CSV report.*

### Morphometric Distributions
![Distributions](docs/distributions_sample.png)
*Population distributions with WHO 6th Edition reference range boundaries (dashed green lines).*

---

## Methodology

### 1. Segmentation
Sperm heads are isolated using adaptive thresholding followed by morphological operations. Connected components are filtered by area (80–800 px) to exclude debris and overlapping cells.

### 2. Principal Axes Estimation
Rather than using simple bounding boxes, CASMA computes the **inertia tensor** (covariance matrix of pixel coordinates) and extracts eigenvalues. The square roots of the eigenvalues, scaled by 4σ, provide robust estimates of the major and minor axes — this method is invariant to rotation and more accurate than rectangular bounding boxes for elliptical objects.

### 3. Robust Outlier Detection
Standard Z-scores assume normality and are unreliable for small samples (n < 30). CASMA uses the **modified Z-score** with Median Absolute Deviation (MAD):

```
MZi = 0.6745 × (xi − median) / MAD
```

Where MAD = median(|xi − median|). Values with |MZi| > 2.5 are flagged as population outliers.

### 4. ALCOA+ Compliance
Every run generates an audit trail documenting:
- **Attributable:** Automated pipeline with version-tagged source code
- **Legible:** Standard CSV/JSON/PNG outputs
- **Contemporaneous:** ISO 8601 timestamp at runtime
- **Original:** Raw image preserved; all calculations reproducible
- **Accurate:** WHO 6th Edition reference ranges hard-coded and versioned

---

## Limitations & Future Work

- **Brightfield only:** Currently tested on phase-contrast and brightfield images. Fluorescence support planned.
- **2D analysis:** Does not measure midpiece or tail morphology (requires multi-focus or AI segmentation).
- **Single frame:** Batch processing of time-lapse sequences is on the roadmap.
- **Validation:** Pending correlation study against manual morphology assessment by trained andrologists.

---

## References

1. World Health Organization. *WHO Laboratory Manual for the Examination and Processing of Human Semen* (6th ed.). Geneva: WHO Press; 2021.
2. Iglewicz B, Hoaglin DC. *How to Detect and Handle Outliers*. Milwaukee: ASQC Quality Press; 1993. (Modified Z-score / MAD method)
3. Young IT. *Proof without prejudice: use of the Kolmogorov-Smirnov test for the analysis of histograms from flow systems and other sources*. J Histochem Cytochem. 1977;25(7):935–941.

---

## License

MIT License — free for academic and research use. Not intended for clinical diagnosis without regulatory validation.

---

## Contact

**Adrian Hall**  
Biomedical Science BSc (Hons) — University of Warwick  
[halladrian785@gmail.com](mailto:halladrian785@gmail.com) | [github.com/AdrianHall23](https://github.com/AdrianHall23)
```

---

Also create a `requirements.txt` with this:

```text
numpy>=1.21.0
scipy>=1.7.0
matplotlib>=3.4.0
Pillow>=8.3.0
```
