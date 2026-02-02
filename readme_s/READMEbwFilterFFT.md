  # Polarization Vector Map bwFilterFFT.py Toolkit Setup


- **NOTE**: As of 1/21/2026 TEMUL Toolkit has not been updated since 2022 and does NOT work with modern versions of certain packages.
This issue may get worse as time goes on, but the required package versions are noted below.

It is highly recommended that a separate Python environment is created to prevent these deprecated package versions to not interfere
with other projects. This can be done through Anaconda or Miniconda3.

## Requirements

- **Python version**: 3.11

- **Scipy version**: 1.10.0
- **Atomap version**: 0.3.4
- **Hyperspy version**: 1.7.5
- **Temul version**: 0.1.7


Download packages beginning with hyperspy -> atomap -> scipy -> temul.

## Setup Instructions

### 1. Add the Input Image

- Locate the image file from downloads:  
  `IFFT of Non-regid alignment 002-2_bw copy.png`
- Copy this file into your PyCharm Project folder.
- In PyCharm, right-click on the image file and select **Copy Path**.
- Open `bwFilterFFT.py`.
- Go to **line 101**, delete the existing file path, and **paste your copied path** in its place.
- Repeat these steps for any image.

### 2. Install Required Packages

#### a. From Dropbox (Manual Import)

Make sure the following Python packages are imported **from the Dropbox folder** provided:

- `Temul`
- `ase`
- `atomap`

> Copy these folders into your PyCharm project directory within the library root.

#### b. From PyCharm (Install via PyCharm UI or pip)

All other dependencies can be installed via PyCharm:

- Go to **PyCharm > Preferences > Project: [Your Project Name] > Python Interpreter**
- Use the `+` button to install required packages (e.g., `numpy`, `scipy`, `matplotlib`, etc.)

---
