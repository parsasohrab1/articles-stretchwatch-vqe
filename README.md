# articles-stretchwatch-vqe
```markdown
# aeticles-stretchwatch-vqe.md

## SRS — استرچ‌واچ: کیفیت نقطه‌ی شروع VQE

### 1. Introduction
- **Purpose:** تعریف ابزاری برای تحلیل کمی کیفیت حالت Hartree-Fock به‌عنوان نقطه‌ی شروع VQE در طول منحنی کشش پیوند مولکول‌های کوچک.
- **Scope:** محاسبه‌ی هامیلتونی کیوبیتی برای H2، HeH+، LiH در فاصله‌های پیوند مختلف، محاسبه‌ی همپوشانی HF با حالت پایه و آنتروپی درهم‌تنیدگی.
- **Definitions:** VQE، Hartree-Fock، bond stretching، overlap، entanglement entropy.

### 2. Overall Description
- **Product Functions:** تولید هامیلتونی در فاصله‌های مختلف، حل دقیق (full CI) برای مرجع، محاسبه‌ی همپوشانی HF و آنتروپی، رسم منحنی‌ها.
- **User Characteristics:** پژوهشگر شیمی کوانتومی با دانش پایتون و PennyLane.
- **Constraints:** محاسبات کلاسیک دقیق برای مولکول‌های کوچک؛ بدون سخت‌افزار کوانتومی.
- **Assumptions:** هامیلتونی‌ها با پایه‌ی STO-3G و تبدیل Jordan-Wigner ساخته می‌شوند.

### 3. Specific Requirements
- **Functional:**
  - FR-01: تعریف هندسه‌ی مولکول و تولید هامیلتونی با PennyLane/PySCF.
  - FR-02: حل دقیق حالت پایه (full configuration interaction) برای مرجع.
  - FR-03: ساخت حالت HF و محاسبه‌ی همپوشانی با حالت پایه.
  - FR-04: محاسبه‌ی آنتروپی درهم‌تنیدگی اوربیتال‌های طبیعی.
  - FR-05: تکرار برای آرایه‌ای از فاصله‌های پیوند و ذخیره‌ی نتایج در CSV.
  - FR-06: تعیین آستانه‌ی فاصله‌ای که پس از آن HF نقطه‌ی شروع نامناسبی می‌شود.
- **External Interfaces:** ورودی JSON با مشخصات مولکول و آرایه‌ی فاصله‌ها؛ خروجی CSV.
- **Performance:** محاسبه برای H2 در ۲۰ فاصله در کمتر از ۱ دقیقه.
- **Design Constraints:** پایتون، PennyLane، PySCF (اختیاری).

### 4. قابلیت ثبت اختراع و نوآوری پیشنهادی
- **ارزیابی اولیه:** قابلیت ثبت اختراع **پایین**. Hartree-Fock به‌عنوان نقطه‌ی شروع VQE و رفتار آن در کشش پیوند در ادبیات و پتنت‌ها به‌طور گسترده پوشش داده شده است.
- **نوآوری پیشنهادی برای SRS:** افزودن **«معیار کیفیت نقطه‌ی شروع»** به‌صورت یک شاخص ترکیبی از همپوشانی و آنتروپی، و استفاده از آن برای **پیشنهاد خودکار سوئیچ به یک حالت مرجع جایگزین** (مثلاً CISD) در فاصله‌های بزرگ. این یک مطالعه‌ی روش‌شناختی است و ادعای پتنت آن بسیار دشوار خواهد بود.
- **اقدام لازم:** SRS را به‌عنوان یک **ابزار تحلیل و بنچمارک** ارائه دهید، نه یک محصول قابل ثبت.

---

## پیوست: کد تولید داده مصنوعی

```python
#!/usr/bin/env python3
"""تولید منحنی انرژی مصنوعی برای H2، HeH+، LiH (مدل مورس تقریبی)."""
import numpy as np
import pandas as pd
import os

def morse_energy(r, De, a, re):
    """انرژی پتانسیل مورس."""
    return De * (1 - np.exp(-a * (r - re))) ** 2 - De

def generate_curve(molecule, r_min=0.3, r_max=4.0, n=40, seed=42):
    rng = np.random.default_rng(seed)

    # پارامترهای تقریبی هر مولکول (صرفاً برای داده مصنوعی)
    params = {
        "H2":   {"De": 4.75, "a": 1.94, "re": 0.74},
        "HeH+": {"De": 2.50, "a": 1.50, "re": 0.77},
        "LiH":  {"De": 2.50, "a": 0.80, "re": 1.60},
    }[molecule]

    r = np.linspace(r_min, r_max, n)
    E_exact = morse_energy(r, params["De"], params["a"], params["re"])
    E_exact += rng.normal(0, 1e-4, n)

    # انرژی HF: در فاصله‌های کوچک نزدیک به دقیق، در فاصله‌های بزرگ واگرا
    delta = r - params["re"]
    E_hf = E_exact + 0.5 * delta**2 * np.exp(-delta) + rng.normal(0, 2e-4, n)

    # همپوشانی HF با حالت پایه (کاهش با کشش)
    overlap = np.exp(-0.3 * np.abs(delta) ** 1.5) + rng.normal(0, 0.01, n)
    overlap = np.clip(overlap, 0, 1)

    # آنتروپی درهم‌تنیدگی (افزایش با کشش)
    entropy = 0.1 + 0.3 * (1 - np.exp(-0.5 * np.abs(delta))) + rng.normal(0, 0.01, n)

    return pd.DataFrame({
        "r": r,
        "E_exact": E_exact,
        "E_HF": E_hf,
        "overlap": overlap,
        "entropy": entropy,
    })

if __name__ == "__main__":
    out_dir = "data/stretchwatch"
    os.makedirs(out_dir, exist_ok=True)
    for mol in ["H2", "HeH+", "LiH"]:
        df = generate_curve(mol, seed=hash(mol) % 1000)
        df.to_csv(f"{out_dir}/{mol}_curve.csv", index=False)
    print(f"منحنی‌های کشش پیوند در {out_dir} تولید شد.")
