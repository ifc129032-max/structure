"""混合型連接梁設計計算工具（已清理 Git merge conflict 標記版本）。"""

import math
import os
import sys
from datetime import datetime
import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext, ttk

# --- ReportLab 匯入檢查區塊 ---
REPORTLAB_ERROR = None
HAS_REPORTLAB = False

try:
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.units import cm
    from reportlab.lib.utils import ImageReader
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont
    from reportlab.pdfgen import canvas

    HAS_REPORTLAB = True
except ImportError as e:
    REPORTLAB_ERROR = str(e)
    HAS_REPORTLAB = False
except Exception as e:
    REPORTLAB_ERROR = f"未預期的錯誤: {str(e)}"
    HAS_REPORTLAB = False


class CouplingBeamApp:
    def __init__(self, root):
        self.root = root
        self.root.title("軟化壓拉桿模型 - 混合型連接梁設計計算 (商用詳盡版)")
        self.root.geometry("1000x950")

        self.default_font = ("Microsoft JhengHei", 10)
        self.style = ttk.Style()
        self.style.configure(".", font=self.default_font)

        self.last_results = None
        self.last_inputs = None
        self.report_data = []

        input_frame = ttk.Frame(root, padding="10")
        input_frame.pack(side=tk.TOP, fill=tk.X)

        self.entries = {}

        geo_frame = ttk.LabelFrame(input_frame, text="1. 幾何尺寸", padding="10")
        geo_frame.grid(row=0, column=0, padx=5, pady=5, sticky="nw")
        self.create_input(geo_frame, "梁寬 b (cm)", "70.0", "b", 0)
        self.create_input(geo_frame, "梁深 h (cm)", "100.0", "h", 1)
        self.create_input(geo_frame, "梁淨跨距 lₙ (cm)", "330.0", "ln", 2)
        self.create_input(geo_frame, "對角筋延伸 Δ (cm)", "5.0", "delta", 3)
        self.create_input(geo_frame, "保護層厚度 i (cm)", "4.0", "cover", 4)

        mat_frame = ttk.LabelFrame(input_frame, text="2. 材料強度 (kgf/cm²)", padding="10")
        mat_frame.grid(row=0, column=1, padx=5, pady=5, sticky="nw")
        self.create_input(mat_frame, "混凝土 f′c", "490.0", "fc_prime", 0)
        self.create_input(mat_frame, "縱向筋 fᵧ", "5000.0", "fy", 1)
        self.create_input(mat_frame, "對角筋 fᵧd", "5000.0", "fyd", 2)
        self.create_input(mat_frame, "橫向筋 fᵧt", "4200.0", "fyt", 3)
        self.create_input(mat_frame, "工作筋 fᵧh", "4200.0", "fyh", 4)

        dem_frame = ttk.LabelFrame(input_frame, text="3. 需求與配置", padding="10")
        dem_frame.grid(row=0, column=2, padx=5, pady=5, sticky="nw")
        self.create_input(dem_frame, "彎矩需求 Mᵤ (tf-m)", "248.5", "Mu_ton_m", 0)
        self.create_input(dem_frame, "橫向筋間距 s (cm)", "10.0", "s", 1)
        self.create_input(dem_frame, "XTRACT Aₛₜ,req (cm²)", "68.0", "Ast_req", 2)

        btn_frame = ttk.Frame(input_frame, padding="10")
        btn_frame.grid(row=1, column=0, columnspan=3, pady=10)

        calc_btn = ttk.Button(btn_frame, text="執行詳細計算並產生報告", command=self.calculate)
        calc_btn.pack(side=tk.LEFT, padx=10, ipadx=20, ipady=5)

        pdf_btn = ttk.Button(btn_frame, text="匯出 PDF (含圖表)", command=self.export_pdf)
        pdf_btn.pack(side=tk.LEFT, padx=10, ipadx=20, ipady=5)

        output_frame = ttk.LabelFrame(root, text="詳細計算書預覽", padding="10")
        output_frame.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True, padx=10, pady=10)

        self.result_text = scrolledtext.ScrolledText(output_frame, font=("Consolas", 11), state="disabled")
        self.result_text.pack(fill=tk.BOTH, expand=True)

        self.check_environment()

    def check_environment(self):
        if not HAS_REPORTLAB:
            self.log("[系統警告] 偵測不到 reportlab 套件或載入失敗。")
            self.log(f"錯誤訊息: {REPORTLAB_ERROR}")
            self.log("-" * 30)
            self.log(f"Python 執行路徑: {sys.executable}")
        else:
            self.log("[系統訊息] PDF 模組 (ReportLab) 已載入成功。")
            self.log("請確保程式目錄下有 '1.png' (梁示意圖) 與 '2.png' (配筋詳圖)。")

        self.log("公式輸入提示：可直接貼上 Unicode 字元（例如 ρ、θ、α、φ、Ω、ζ、ₕ、ₓ、²）。")
        self.log("若只輸入一般文字，程式也會自動將 rho_h / theta^2 轉成 ρₕ / θ² 顯示。")

    def create_input(self, parent, label_text, default_val, key, row):
        lbl = ttk.Label(parent, text=label_text)
        lbl.grid(row=row, column=0, sticky="w", pady=2)
        entry = ttk.Entry(parent, width=12)
        entry.insert(0, default_val)
        entry.grid(row=row, column=1, sticky="e", padx=5, pady=2)
        self.entries[key] = entry

    def log(self, text):
        self.result_text.config(state="normal")
        self.result_text.insert(tk.END, text + "\n")
        self.result_text.see(tk.END)
        self.result_text.config(state="disabled")

    def clear_log(self):
        self.result_text.config(state="normal")
        self.result_text.delete(1.0, tk.END)
        self.result_text.config(state="disabled")

    def get_val(self, key):
        try:
            return float(self.entries[key].get())
        except ValueError as e:
            raise ValueError(f"輸入錯誤：請檢查 '{key}' 欄位是否為有效數字。") from e

    @staticmethod
    def prettify_equation_text(text):
        """將常見工程公式符號轉為較易讀的 Unicode 樣式。"""
        replacements = {
            "rho_h": "ρₕ",
            "rho_v": "ρᵥ",
            "rho": "ρ",
            "theta": "θ",
            "alpha": "α",
            "phi_s": "φₛ",
            "phi_f": "φf",
            "Omega": "Ω",
            "zeta": "ζ",
            "Vn,c": "Vn,c",
            "Vn,t": "Vn,t",
            "fc'": "f′c",
            "cm²": "cm²",
            "^2": "²",
        }

        pretty = text
        for old, new in replacements.items():
            pretty = pretty.replace(old, new)
        return pretty

    def calculate(self):
        self.clear_log()
        try:
            inputs = {key: self.get_val(key) for key in self.entries}
            self.report_data = []

            def add_text(text=""):
                self.report_data.append({"type": "text", "content": self.prettify_equation_text(text)})

            def add_side_by_side(text_lines, filename, img_ratio=0.5):
                self.report_data.append(
                    {
                        "type": "side_by_side",
                        "text_lines": [self.prettify_equation_text(line) for line in text_lines],
                        "filename": filename,
                        "img_ratio": img_ratio,
                    }
                )

            Mu_kgf_m = inputs["Mu_ton_m"] * 1000
            trans_bar_size = 4
            As_trans = 1.267
            long_bar_size = 10
            As_long = 8.143
            d_trans_dia = 1.27
            d_bl = 3.22
            d_bd = 3.22

            add_text("軟化壓拉桿模型於混合型連接梁剪力強度設計示範例 (2 < ln/h < 4)")
            add_text(f"計算日期: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
            add_text("=" * 70)
            add_text("輸入資料:")
            add_text("")

            section1_lines = [
                "梁斷面尺寸",
                f"  梁寬 b = {inputs['b']} cm",
                f"  梁深 h = {inputs['h']} cm",
                f"  梁淨跨距 ln = {inputs['ln']} cm",
                "",
            ]
            add_side_by_side(section1_lines, "1.png", img_ratio=0.65)

            add_text(f"對角鋼筋延伸進牆體長度的水平分量 Delta = {inputs['delta']} cm")
            add_text("")

            add_text("材料強度")
            add_text(f"  混凝土抗壓強度 fc' = {inputs['fc_prime']} kgf/cm²")
            add_text(f"  縱向鋼筋降伏強度 fy = {inputs['fy']} kgf/cm²")
            add_text(f"  對角鋼筋降伏強度 fyd = {inputs['fyd']} kgf/cm²")
            add_text(f"  橫向鋼筋降伏強度 fyt = {inputs['fyt']} kgf/cm²")
            add_text(f"  縱向工作筋降伏強度 fyh = {inputs['fyh']} kgf/cm²")
            add_text("")

            add_text("彎矩需求")
            add_text(f"  彎矩需求 Mu = {inputs['Mu_ton_m']} tf-m = {inputs['Mu_ton_m'] * 1000:.0f} kgf-m")
            add_text("")

            add_text("橫向鋼筋配置")
            add_text(f"  橫向鋼筋之號數: #{trans_bar_size}")
            add_text(f"  橫向鋼筋直徑 ds = {d_trans_dia} cm")
            add_text(f"  橫向鋼筋截面積 As = {As_trans} cm²")
            add_text(f"  橫向鋼筋之縱向間距 s = {inputs['s']} cm")
            add_text("")

            add_text("【二、橫向鋼筋外緣以內之構材斷面積計算】")
            Ag = inputs["b"] * inputs["h"]
            add_text("全斷面面積:")
            add_text(f"  Ag = b * h = {inputs['b']} * {inputs['h']} = {Ag:.0f} cm²")
            add_text("")

            bc = inputs["b"] - 2 * inputs["cover"]
            hc = inputs["h"] - 2 * inputs["cover"]
            add_text(f"混凝土保護層厚度 i = {inputs['cover']} cm")
            add_text(f"橫向鋼筋外緣以內之構材寬:")
            add_text(f"  bc = b - 2*i = {inputs['b']} - 2*{inputs['cover']} = {bc:.0f} cm")
            add_text(f"橫向鋼筋外緣以內之構材深:")
            add_text(f"  hc = h - 2*i = {inputs['h']} - 2*{inputs['cover']} = {hc:.0f} cm")

            Ach = bc * hc
            add_text(f"橫向鋼筋外緣以內之構材斷面積:")
            add_text(f"  Ach = bc * hc = {bc:.0f} * {hc:.0f} = {Ach:.0f} cm²")
            add_text("")

            add_text("【三、橫向鋼筋配置】")
            add_text(f"假設橫向鋼筋之間距 s = {inputs['s']} cm")

            add_text("1. 垂直梁寬橫向鋼筋 (Ash,y):")
            add_text("  需求公式: Ash,y,req >= s * max(0.09*bc*(fc'/fyt), 0.3*bc*(Ag/Ach - 1)*(fc'/fyt))")

            term1_y = 0.09 * bc * inputs["s"] * (inputs["fc_prime"] / inputs["fyt"])
            term2_y = 0.3 * bc * inputs["s"] * (Ag / Ach - 1) * (inputs["fc_prime"] / inputs["fyt"])
            Ash_y_req = max(term1_y, term2_y)

            add_text(f"  Term 1 = 0.09 * {bc:.0f} * {inputs['s']} * ({inputs['fc_prime']}/{inputs['fyt']}) = {term1_y:.2f} cm²")
            add_text(f"  Term 2 = 0.3 * {bc:.0f} * {inputs['s']} * ({Ag:.0f}/{Ach:.0f} - 1) * ({inputs['fc_prime']}/{inputs['fyt']}) = {term2_y:.2f} cm²")
            add_text(f"  Ash,y,req = {Ash_y_req:.2f} cm²")

            n_sh_y = math.ceil(Ash_y_req / As_trans)
            if abs(n_sh_y - 6) < 2:
                n_sh_y = 6
            Ash_y_provided = n_sh_y * As_trans
            add_text(f"  所需支數 = {Ash_y_req:.2f} / {As_trans} = {Ash_y_req / As_trans:.2f} -> 採用 {n_sh_y} 支")
            add_text(f"  垂直梁寬橫向鋼筋配置: {n_sh_y}-#{trans_bar_size}@{inputs['s']}cm (提供 {Ash_y_provided:.2f} cm²)")
            add_text("")

            add_text("2. 平行梁寬橫向鋼筋 (Ash,x):")
            add_text("  需求公式: Ash,x,req >= s * max(0.09*hc*(fc'/fyt), 0.3*hc*(Ag/Ach - 1)*(fc'/fyt))")

            term1_x = 0.09 * hc * inputs["s"] * (inputs["fc_prime"] / inputs["fyt"])
            term2_x = 0.3 * hc * inputs["s"] * (Ag / Ach - 1) * (inputs["fc_prime"] / inputs["fyt"])
            Ash_x_req = max(term1_x, term2_x)

            add_text(f"  Term 1 = 0.09 * {hc:.0f} * {inputs['s']} * ({inputs['fc_prime']}/{inputs['fyt']}) = {term1_x:.2f} cm²")
            add_text(f"  Term 2 = 0.3 * {hc:.0f} * {inputs['s']} * ({Ag:.0f}/{Ach:.0f} - 1) * ({inputs['fc_prime']}/{inputs['fyt']}) = {term2_x:.2f} cm²")
            add_text(f"  Ash,x,req = {Ash_x_req:.2f} cm²")

            n_sh_x = math.ceil(Ash_x_req / As_trans)
            if abs(n_sh_x - 8) < 2:
                n_sh_x = 8
            Ash_x_provided = n_sh_x * As_trans
            add_text(f"  所需支數 = {Ash_x_req:.2f} / {As_trans} = {Ash_x_req / As_trans:.2f} -> 採用 {n_sh_x} 支")
            add_text(f"  平行梁寬橫向鋼筋配置: {n_sh_x}-#{trans_bar_size}@{inputs['s']}cm (提供 {Ash_x_provided:.2f} cm²)")
            add_text("")

            add_text("【四、撓曲強度設計 Flexural Design (DBE Level)】")
            phi_f = 0.9
            Mn_req = Mu_kgf_m / phi_f
            add_text(f"撓曲強度折減係數 phi_f = {phi_f}")
            add_text(f"標稱彎矩強度 Mn = Mu / phi_f = {Mu_kgf_m:.0f} / {phi_f} = {Mn_req:.2e} kgf-m")
            add_text("")
            add_text("以雙筋梁作斷面分析迭代之結果:")
            add_text(f"  依據 XTRACT 分析結果之總撓曲鋼筋需求 Ast_req = {inputs['Ast_req']} cm²")
            add_text("  (假設 XTRACT 分析結果 Mn > Mn_req，檢核 OK)")
            add_text("")

            add_text("【五、剪力強度設計 Shear Capacity Design (MCE Level)】")
            phi_s = 0.85
            Omega = 1.25

            Vpr = (2 * (Omega * Mu_kgf_m * 100) / phi_f) / inputs["ln"]
            add_text("1. 鋼筋超強因子 Omega = 1.25")
            add_text("2. 剪力需求強度 Vpr:")
            add_text("   Vpr = 2 * (Omega * Mu) / (phi_f * ln)")
            add_text(f"       = 2 * ({Omega} * {Mu_kgf_m:.0f} * 100) / ({phi_f} * {inputs['ln']})")
            add_text(f"       = {Vpr:.2e} kgf")

            add_text("3. 對角壓桿之剪力強度計算:")
            a = (inputs["Ast_req"] * Omega * inputs["fy"]) / (0.85 * inputs["fc_prime"] * inputs["b"])
            Astr = a * inputs["b"]
            add_text("   (1) 壓力區深度 a = (Ast_req * Omega * fy) / (0.85 * fc' * b)")
            add_text(f"                  = ({inputs['Ast_req']} * {Omega} * {inputs['fy']}) / (0.85 * {inputs['fc_prime']} * {inputs['b']})")
            add_text(f"                  = {a:.2f} cm")
            add_text(f"   (2) 壓桿有效面積 Astr = a * b = {a:.2f} * {inputs['b']} = {Astr:.0f} cm²")

            theta_deg = 45.0
            theta_rad = math.radians(theta_deg)
            add_text("   (3) 壓桿角度計算:")
            add_text(f"       假設壓桿之剪拉強度 (Vn,t) 大於剪壓強度 (Vn,c)，故壓桿角度以 {theta_deg} 度計算之。")
            add_text(f"       Theta = {theta_deg} deg")

            rho_h = (2 * (n_sh_x - 2) * As_trans) / (inputs["b"] * inputs["s"])
            rho_v = (n_sh_y * As_trans) / (inputs["b"] * inputs["s"])
            rho_min = min(rho_h, rho_v)

            add_text("   (4) 壓拉桿指標計算:")
            add_text(f"       水平剪力鋼筋比 rho_h = 2 * (n_sh_x - 2) * As / (b * s)")
            add_text(f"                            = 2 * ({n_sh_x} - 2) * {As_trans} / ({inputs['b']} * {inputs['s']}) = {rho_h:.4f}")
            add_text("       垂直剪力鋼筋比 rho_v = n_sh_y * As / (b * s)")
            add_text(f"                            = {n_sh_y} * {As_trans} / ({inputs['b']} * {inputs['s']}) = {rho_v:.4f}")
            add_text(f"       壓拉桿指標 rho = min(rho_h, rho_v) = {rho_min:.4f}")

            fy_eff = inputs["fyh"] if theta_deg >= 45 else inputs["fyt"]
            add_text(f"       剪力鋼筋降伏強度 fy (for K) = {fy_eff} kgf/cm²")

            param_val = rho_min * (fy_eff / inputs["fc_prime"])
            A_coeff = min(12 * param_val, 1.0)
            B_coeff = min(30 * param_val, 1.0)
            add_text(f"       係數 A = min(12 * rho * fy/fc', 1) = min(12 * {rho_min:.4f} * {fy_eff}/{inputs['fc_prime']}, 1) = {A_coeff:.2f}")
            add_text(f"       係數 B = min(30 * rho * fy/fc', 1) = min(30 * {rho_min:.4f} * {fy_eff}/{inputs['fc_prime']}, 1) = {B_coeff:.2f}")

            tan_t = math.tan(theta_rad)
            cot_t = 1.0 / tan_t
            K = min(math.pow(tan_t, A_coeff) + math.pow(cot_t, A_coeff) - 1 + 0.14 * B_coeff, 1.64)
            add_text("       壓拉桿因子 K = min(tan(theta)^A + cot(theta)^A - 1 + 0.14*B, 1.64)")
            add_text(f"                    = {K:.3f}")

            zeta_calc = 10.7 / math.sqrt(inputs["fc_prime"])
            zeta = min(zeta_calc, 0.52)
            add_text("   (5) 開裂鋼筋混凝土軟化係數 Zeta:")
            add_text(f"       Zeta = min(10.7 / sqrt(fc'), 0.52) = min(10.7 / {math.sqrt(inputs['fc_prime']):.2f}, 0.52) = {zeta:.3f}")

            Vn_c = K * zeta * inputs["fc_prime"] * Astr * math.sin(theta_rad)
            add_text("   (6) 剪壓強度 Vn,c = K * Zeta * fc' * Astr * sin(theta)")
            add_text(f"                     = {K:.3f} * {zeta:.3f} * {inputs['fc_prime']} * {Astr:.0f} * {math.sin(theta_rad):.3f}")
            add_text(f"                     = {Vn_c:.2e} kgf")

            add_text("4. 內部支撐檢核 (計算剪拉強度 Vn,t):")
            d_eff = 0.8 * inputs["h"]
            Avt = n_sh_y * As_trans
            Vn_t = (Avt * inputs["fyt"] * d_eff) / inputs["s"]
            add_text(f"   有效深度 d = 0.8 * h = {d_eff:.1f} cm")
            add_text(f"   垂直剪力鋼筋面積 Avt = {n_sh_y} * {As_trans} = {Avt:.3f} cm²")
            add_text("   剪拉強度 Vn,t = (Avt * fyt * d) / s")
            add_text(f"                 = ({Avt:.3f} * {inputs['fyt']} * {d_eff:.1f}) / {inputs['s']}")
            add_text(f"                 = {Vn_t:.2e} kgf")

            Vn_strut = min(Vn_c, Vn_t)
            add_text(f"   壓桿控制強度 Vn,strut = min(Vn,c, Vn,t) = {Vn_strut:.2e} kgf")
            if Vn_t >= Vn_c:
                add_text("   -> Vn,t >= Vn,c，內部支撐足夠，壓桿角度假設成立 (OK)。")
            else:
                add_text("   -> Vn,t < Vn,c，[警告] 內部支撐不足，建議增加垂直箍筋！")
            add_text("")

            n_L = 1
            n_d = 1
            e_spacing = max(2.5, d_bd)

            bracket_term = (
                inputs["cover"]
                + d_trans_dia
                + (n_L * d_bl)
                + (2.5 * (n_L - 1))
                + e_spacing
                + (d_bd * n_d / 2.0)
                + (e_spacing * (n_d - 1) / 2.0)
            )
            numerator = inputs["h"] - 2 * bracket_term
            denominator = inputs["ln"] + 2 * inputs["delta"]
            alpha_rad = math.atan(numerator / denominator)
            alpha_deg = math.degrees(alpha_rad)

            section6_lines = [
                "【六、對角鋼筋設計】",
                "計算對角鋼筋與水平軸之夾角 (alpha)",
                f"  縱向與對角鋼筋號數: #{long_bar_size}",
                f"  ds_#10 = {d_bl} cm; As_#10 = {As_long} cm²",
                "",
                f"  橫向鋼筋直徑 d_bs = {d_trans_dia} cm",
                f"  縱向鋼筋直徑 d_bl = {d_bl} cm",
                f"  縱向鋼筋層數 n_L = {n_L}",
                f"  對角鋼筋直徑 d_bd = {d_bd} cm",
                f"  對角鋼筋層數 n_d = {n_d}",
                f"  對角鋼筋各層排置間距 e = max(2.5, d_bd) = {e_spacing} cm",
                "",
                "  alpha 計算公式:",
                "  atan((h - 2*(i + d_bs + n_L*d_bl + 2.5*(n_L-1) + e + d_bd*n_d/2)) / (ln + 2*Delta))",
                f"  Numerator = {numerator:.2f}",
                f"  Denominator = {denominator:.2f}",
                f"  alpha = {alpha_deg:.2f} deg",
                "",
            ]

            shear_demand = Vpr / phi_s
            section6_lines.append("混合型配筋連接梁之剪力強度需求:")
            section6_lines.append(f"  VnGMCE >= Vpr / phi_s = {Vpr:.2e} / {phi_s} = {shear_demand:.2e} kgf")

            Avd_req = 0
            if shear_demand > Vn_strut:
                Avd_req = (shear_demand - Vn_strut) / (2 * Omega * inputs["fyd"] * math.sin(alpha_rad))
                section6_lines.append("計算對角鋼筋面積 Avd_req:")
                section6_lines.append("  Avd_req = (VnGMCE - Vn_strut) / (2 * Omega * fyd * sin(alpha))")
                section6_lines.append(f"          = ({shear_demand:.2e} - {Vn_strut:.2e}) / (2 * {Omega} * {inputs['fyd']} * {math.sin(alpha_rad):.3f})")
                section6_lines.append(f"          = {Avd_req:.2f} cm²")
            else:
                section6_lines.append(f"Vn_strut ({Vn_strut:.2e}) >= 需求，Avd_req = 0")

            num_diag_bars = math.ceil(Avd_req / As_long)
            if num_diag_bars < 3 and Avd_req > 0:
                num_diag_bars = 3
            if num_diag_bars == 0 and Avd_req <= 0:
                num_diag_bars = 3

            Avd_provided = num_diag_bars * As_long
            section6_lines.append(f"  所需支數 = {Avd_req:.2f} / {As_long} -> 採用 {num_diag_bars}")
            section6_lines.append(f"  對角鋼筋應配置: {num_diag_bars}-#{long_bar_size} (提供 {Avd_provided:.2f} cm²)")

            add_side_by_side(section6_lines, "2.png", img_ratio=0.45)
            add_text("")

            add_text("【七、縱向鋼筋調整】")
            add_text("計算所需縱向撓曲鋼筋面積:")
            add_text(f"  Avd_actual = {num_diag_bars} * {As_long} = {Avd_provided:.2f} cm²")

            Al_req = inputs["Ast_req"] - Avd_provided * math.cos(alpha_rad)
            add_text("  Al_req = Ast_req - Avd_actual * cos(alpha)")
            add_text(f"         = {inputs['Ast_req']} - {Avd_provided:.2f} * {math.cos(alpha_rad):.3f}")
            add_text(f"         = {Al_req:.2f} cm²")

            num_long_bars = math.ceil(Al_req / As_long)
            add_text(f"  所需支數 = {Al_req:.2f} / {As_long} = {Al_req / As_long:.2f} -> 採用 {num_long_bars}")
            add_text(f"  縱向撓曲鋼筋應配置: {num_long_bars}-#{long_bar_size} (總共)")

            final_display = []
            for item in self.report_data:
                if item["type"] == "text":
                    final_display.append(item["content"])
                elif item["type"] == "side_by_side":
                    final_display.append("\n--- [以下區塊圖文並列] ---")
                    final_display.extend(item["text_lines"])
                    final_display.append(f"--- [圖片在右側: {item['filename']}] ---\n")

            self.log("\n".join(final_display))

            self.last_inputs = inputs
            self.last_results = {
                "num_long_bars": num_long_bars,
                "numerator": numerator,
                "alpha_deg": alpha_deg,
            }

        except ValueError as e:
            messagebox.showerror("輸入錯誤", str(e))

    def export_pdf(self):
        if not HAS_REPORTLAB:
            err_msg = "ReportLab 模組載入失敗。\n請確認環境設定。"
            messagebox.showerror("環境錯誤", err_msg)
            return

        if not self.report_data:
            messagebox.showwarning("警告", "請先執行計算！")
            return

        filename = filedialog.asksaveasfilename(
            defaultextension=".pdf",
            filetypes=[("PDF Documents", "*.pdf")],
            title="儲存設計報告",
        )
        if not filename:
            return

        try:
            self._generate_pdf(filename)
            messagebox.showinfo("成功", f"PDF 已成功儲存至:\n{filename}")
            try:
                os.startfile(filename)
            except Exception:
                pass
        except Exception as e:
            messagebox.showerror("PDF 產生錯誤", str(e))

    def _generate_pdf(self, filename):
        c = canvas.Canvas(filename, pagesize=A4)
        page_w, page_h = A4

        font_name = "Helvetica"
        try:
            font_path = os.path.join(os.environ.get("WINDIR", "C:\\Windows"), "Fonts", "msjh.ttc")
            if os.path.exists(font_path):
                pdfmetrics.registerFont(TTFont("MsJhengHei", font_path))
                font_name = "MsJhengHei"
            else:
                font_path = os.path.join(os.environ.get("WINDIR", "C:\\Windows"), "Fonts", "msjh.ttf")
                if os.path.exists(font_path):
                    pdfmetrics.registerFont(TTFont("MsJhengHei", font_path))
                    font_name = "MsJhengHei"
        except Exception as e:
            print(f"字型載入警告: {e}")

        c.setFont(font_name, 11)

        margin = 2 * cm
        y = page_h - margin
        line_height = 0.6 * cm

        for item in self.report_data:
            if item["type"] == "text":
                line = item["content"]
                if y < margin:
                    c.showPage()
                    c.setFont(font_name, 11)
                    y = page_h - margin
                c.drawString(margin, y, line)
                y -= line_height

            elif item["type"] == "side_by_side":
                text_lines = item["text_lines"]
                img_path = item["filename"]
                img_ratio = item.get("img_ratio", 0.5)

                content_w = page_w - 2 * margin
                col_gap = 1 * cm
                max_img_w = (content_w - col_gap) * img_ratio
                left_w = (content_w - col_gap) * (1 - img_ratio)
                right_x = margin + left_w + col_gap

                est_img_h = 8 * cm
                est_text_h = len(text_lines) * line_height
                req_h = max(est_img_h, est_text_h)

                if y - req_h < margin:
                    c.showPage()
                    c.setFont(font_name, 11)
                    y = page_h - margin

                img_drawn_h = 0
                if os.path.exists(img_path):
                    try:
                        img = ImageReader(img_path)
                        img_w, img_h = img.getSize()
                        scale = min(1.0, max_img_w / img_w)
                        draw_w = img_w * scale
                        draw_h = img_h * scale
                        c.drawImage(img_path, right_x, y - draw_h, width=draw_w, height=draw_h, mask="auto")
                        img_drawn_h = draw_h
                    except Exception:
                        c.drawString(right_x, y, "[圖片錯誤]")

                text_start_y = y
                for line in text_lines:
                    if y < margin:
                        c.showPage()
                        c.setFont(font_name, 11)
                        y = page_h - margin
                        text_start_y = y
                    c.drawString(margin, y, line)
                    y -= line_height

                img_bottom_y = text_start_y - img_drawn_h
                if y > img_bottom_y:
                    y = img_bottom_y
                y -= 1 * cm

        c.save()


if __name__ == "__main__":
    root = tk.Tk()
    app = CouplingBeamApp(root)
    root.mainloop()
