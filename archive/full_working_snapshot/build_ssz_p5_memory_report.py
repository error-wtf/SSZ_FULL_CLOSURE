from docx import Document
from docx.shared import Inches, Pt, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT, WD_CELL_VERTICAL_ALIGNMENT
from docx.enum.section import WD_SECTION_START
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.enum.style import WD_STYLE_TYPE
from pathlib import Path
import pandas as pd
import numpy as np

B=Path('/mnt/data')
OUT=B/'SSZ_P5_HSVT_FINITE_L_TECHNICAL_MEMORY_2026-09-16.docx'

def shade(cell, fill):
    tcPr=cell._tc.get_or_add_tcPr()
    shd=tcPr.find(qn('w:shd'))
    if shd is None:
        shd=OxmlElement('w:shd'); tcPr.append(shd)
    shd.set(qn('w:fill'), fill)

def set_cell_margins(cell, top=80, start=100, bottom=80, end=100):
    tc=cell._tc; tcPr=tc.get_or_add_tcPr(); tcMar=tcPr.first_child_found_in('w:tcMar')
    if tcMar is None:
        tcMar=OxmlElement('w:tcMar'); tcPr.append(tcMar)
    for m,v in [('top',top),('start',start),('bottom',bottom),('end',end)]:
        node=tcMar.find(qn('w:'+m))
        if node is None:
            node=OxmlElement('w:'+m); tcMar.append(node)
        node.set(qn('w:w'),str(v)); node.set(qn('w:type'),'dxa')

def set_repeat_table_header(row):
    trPr=row._tr.get_or_add_trPr(); tblHeader=OxmlElement('w:tblHeader'); tblHeader.set(qn('w:val'),'true'); trPr.append(tblHeader)

def add_page_number(paragraph):
    paragraph.alignment=WD_ALIGN_PARAGRAPH.RIGHT
    run=paragraph.add_run('Page ')
    fldChar1=OxmlElement('w:fldChar'); fldChar1.set(qn('w:fldCharType'),'begin')
    instrText=OxmlElement('w:instrText'); instrText.set(qn('xml:space'),'preserve'); instrText.text=' PAGE '
    fldChar2=OxmlElement('w:fldChar'); fldChar2.set(qn('w:fldCharType'),'end')
    run._r.append(fldChar1); run._r.append(instrText); run._r.append(fldChar2)

def keep_with_next(p):
    pPr=p._p.get_or_add_pPr(); k=OxmlElement('w:keepNext'); pPr.append(k)

def set_keep_together(p):
    pPr=p._p.get_or_add_pPr(); k=OxmlElement('w:keepLines'); pPr.append(k)

doc=Document()
sec=doc.sections[0]
sec.top_margin=Cm(1.65); sec.bottom_margin=Cm(1.55); sec.left_margin=Cm(1.75); sec.right_margin=Cm(1.75)
sec.header_distance=Cm(0.65); sec.footer_distance=Cm(0.65)

styles=doc.styles
styles['Normal'].font.name='Arial'; styles['Normal'].font.size=Pt(9.4)
styles['Normal'].paragraph_format.space_after=Pt(4.3); styles['Normal'].paragraph_format.line_spacing=1.08
for sname,size,color in [('Title',24,'17365D'),('Heading 1',15,'17365D'),('Heading 2',12,'2F5597'),('Heading 3',10.5,'365F91')]:
    st=styles[sname]; st.font.name='Arial'; st.font.size=Pt(size); st.font.color.rgb=RGBColor.from_string(color); st.font.bold=True
    st.paragraph_format.space_before=Pt(9); st.paragraph_format.space_after=Pt(4)
    if sname!='Title': st.paragraph_format.keep_with_next=True
# custom styles
for name, font, size, color, bold in [
    ('CodeBlock','Courier New',8.2,'202020',False),
    ('Equation','Cambria Math',10.2,'000000',False),
    ('Small','Arial',7.8,'404040',False),
    ('Callout','Arial',9.2,'202020',False),
]:
    if name not in styles:
        st=styles.add_style(name,WD_STYLE_TYPE.PARAGRAPH)
    else: st=styles[name]
    st.font.name=font; st.font.size=Pt(size); st.font.color.rgb=RGBColor.from_string(color); st.font.bold=bold
    st.paragraph_format.space_after=Pt(3)

# header/footer
h=sec.header.paragraphs[0]; h.text='SSZ P5 / Horndeski-SVT finite-l closure - technical memory - 16.09.2026'; h.style='Small'
add_page_number(sec.footer.paragraphs[0]); sec.footer.paragraphs[0].style='Small'

# Helpers

def p(text='', style=None, bold_prefix=None, italic=False):
    par=doc.add_paragraph(style=style)
    if bold_prefix and text.startswith(bold_prefix):
        r=par.add_run(bold_prefix); r.bold=True; par.add_run(text[len(bold_prefix):])
    else:
        r=par.add_run(text); r.italic=italic
    return par

def h1(t): return doc.add_heading(t,level=1)
def h2(t): return doc.add_heading(t,level=2)
def h3(t): return doc.add_heading(t,level=3)
def eq(t):
    par=doc.add_paragraph(style='Equation'); par.alignment=WD_ALIGN_PARAGRAPH.CENTER; par.add_run(t); return par

def code(lines):
    par=doc.add_paragraph(style='CodeBlock'); par.paragraph_format.left_indent=Cm(0.45); par.paragraph_format.right_indent=Cm(0.3)
    pPr=par._p.get_or_add_pPr(); shd=OxmlElement('w:shd'); shd.set(qn('w:fill'),'F2F2F2'); pPr.append(shd)
    par.add_run(lines); return par

def bullet(text, level=0):
    par=doc.add_paragraph(style='List Bullet' if level==0 else 'List Bullet 2'); par.add_run(text); return par

def numbered(text):
    par=doc.add_paragraph(style='List Number'); par.add_run(text); return par

def callout(title, text, fill='EAF2F8'):
    table=doc.add_table(rows=1, cols=1); table.alignment=WD_TABLE_ALIGNMENT.CENTER; table.autofit=True
    cell=table.cell(0,0); shade(cell,fill); set_cell_margins(cell,120,160,120,160)
    q=cell.paragraphs[0]; rr=q.add_run(title); rr.bold=True; rr.font.color.rgb=RGBColor.from_string('17365D')
    q.add_run('\n'+text)
    return table

def table(headers, rows, widths=None, font=7.8):
    tbl=doc.add_table(rows=1, cols=len(headers)); tbl.style='Table Grid'; tbl.alignment=WD_TABLE_ALIGNMENT.CENTER
    hdr=tbl.rows[0]; set_repeat_table_header(hdr)
    for j,x in enumerate(headers):
        hdr.cells[j].text=str(x); shade(hdr.cells[j],'D9EAF7'); hdr.cells[j].vertical_alignment=WD_CELL_VERTICAL_ALIGNMENT.CENTER
        for r in hdr.cells[j].paragraphs[0].runs: r.bold=True; r.font.size=Pt(font)
        set_cell_margins(hdr.cells[j])
    for row in rows:
        cells=tbl.add_row().cells
        for j,x in enumerate(row):
            cells[j].text=str(x); cells[j].vertical_alignment=WD_CELL_VERTICAL_ALIGNMENT.CENTER; set_cell_margins(cells[j])
            for par in cells[j].paragraphs:
                par.paragraph_format.space_after=Pt(0)
                for r in par.runs: r.font.size=Pt(font)
    if widths:
        for row in tbl.rows:
            for j,w in enumerate(widths): row.cells[j].width=Cm(w)
    return tbl

def add_picture(path, width_cm=16.0, caption=None):
    par=doc.add_paragraph(); par.alignment=WD_ALIGN_PARAGRAPH.CENTER
    par.add_run().add_picture(str(path), width=Cm(width_cm))
    if caption:
        cp=doc.add_paragraph(caption, style='Small'); cp.alignment=WD_ALIGN_PARAGRAPH.CENTER; cp.runs[0].italic=True
    return par

# Title page
par=doc.add_paragraph(); par.alignment=WD_ALIGN_PARAGRAPH.CENTER; par.paragraph_format.space_before=Pt(45)
r=par.add_run('SSZ P5 / Horndeski-SVT\nFinite-l Closure and QNM Threshold'); r.bold=True; r.font.size=Pt(24); r.font.color.rgb=RGBColor.from_string('17365D')
par=doc.add_paragraph(); par.alignment=WD_ALIGN_PARAGRAPH.CENTER
r=par.add_run('Detailed technical memory and reproducibility report'); r.bold=True; r.font.size=Pt(15); r.font.color.rgb=RGBColor.from_string('2F5597')
par=doc.add_paragraph(); par.alignment=WD_ALIGN_PARAGRAPH.CENTER; par.paragraph_format.space_before=Pt(18)
par.add_run('Frozen status after the final variable-G4 outer candidate audit\n16 September 2026').font.size=Pt(11)
par=doc.add_paragraph(); par.alignment=WD_ALIGN_PARAGRAPH.CENTER; par.paragraph_format.space_before=Pt(20)
par.add_run('Purpose: external technical memory for later reconstruction by Carmen and ChatGPT.\nThis is not a publication claim; it records what was proved, what was merely constructive, what was superseded, and exactly where the final explicit representative failed.').font.size=Pt(9.5)
callout('STOP CHECKPOINT', 'Do not restart ad-hoc representative hunting from scratch. The 41/41 coefficient language, derivative service, constraint maps, Zhang-Kase emitter, Maxwell-Horndeski emitter, and profile reducer are already closed. The unresolved task is a targeted global finite-l representative-selection/control problem.', 'FFF2CC')
doc.add_page_break()

# Executive summary
h1('1. Executive status - the short version')
p('The most important distinction in the whole project is that several different meanings of "closure" were used at different stages. The apparent contradiction between "41/41 closed" and the final finite-l failure disappears once these levels are separated.')
rows=[
('P5 frozen geometry','PASS','Regular, horizonless, two light rings resolved.'),
('Local sector principal witnesses','PASS','Exterior Horndeski, central genuine-SVT, regular Horndeski core.'),
('Constructive local same-action handovers','PASS','C-infinity existence with background re-solving and transverse jet controls.'),
('Unreduced 41-slot coefficient language','CLOSED','Complete a,b,c,d,e,v tables for selected local/sector representatives.'),
('JET9D8 radial derivative service','PASS','Needed higher radial jets and product rules regress.'),
('Full algebraic even constraint maps','PASS','All common pivots nonzero on certified domains.'),
('Profile-aware K,R,G,S,M reducer','PASS','Includes product-rule/IBP terms and antisymmetric S block.'),
('Zhang-Kase Appendix-A emitter','PASS','Complete 41-slot emitter, regression-controlled.'),
('Variable-G4(phi) Maxwell-Horndeski emitter','PASS','Strong-H 39/41 regression after c6,d3,e4 parser fixes.'),
('Global canonical principal K/R/G stream','PASS as principal witness','45,617-point center-to-infinity characteristic stream.'),
('Variable-G4 outer finite-l candidate','REJECTED','sym(K) is negative over much of the band; moreover the variable-G4 transplant was not followed by a fresh full background re-solve.'),
('Certified global finite-l same-action member','OPEN','No center-to-infinity on-shell member has yet passed the complete finite-l gates.'),
('Coupled QNM spectrum','BLOCKED','Must not be computed/claimed until a certified global finite-l member exists.')]
table(['Object','Status','Meaning'],rows,widths=[5.2,3.0,9.2],font=7.7)
callout('Bottom line', 'The statement "41/41 is closed" remains correct. The Sep-16 variable-G4 outer construction is a rejected diagnostic candidate, not a certified global same-action member: after importing the healthier Horndeski principal carrier, the complete background equations were not re-solved for that modified member. Its negative sym(K) is enough to reject the candidate, but it is not a no-go for the still-open global finite-l same-action problem.', 'E2F0D9')

h2('1.1 Kinetic result of the final variable-G4 candidate')
audit=pd.read_csv(B/'ssz_p5_OUTER_FINAL_SAME_ACTION_OPERATOR_AUDIT_2026-09-16.csv')
rows=[]
for _,r in audit.iterrows():
    rows.append((int(r.L),f'{r.min_eig_K:.6g}',f'{int(r.negative_K_rows)} / 2201',f'{r.min_radial_c2:.6g}',int(r.negative_c2_rows),f'{r.min_abs_Dh1:.6g}',f'{r.min_abs_DeltaV:.6g}',f'{r.min_abs_2Lv9:.6g}'))
table(['L','min eig(K)','negative K rows','min radial c^2','neg c^2','min |Dh1|','min |DeltaV|','min |2Lv9|'],rows,font=7.1)
p('Where K is positive, the radial generalized eigenvalues in this audit remain positive. That does not rescue a ghost-sign failure of K itself. The kinetic gate is logically prior to a physical QNM interpretation.')
add_picture(B/'fig_final_K_eig.png',15.5,'Figure 1. Minimum eigenvalue of sym(K) for the final variable-G4 outer candidate, L=6, 42 and 1000.')
add_picture(B/'fig_final_K_eig_zoom.png',15.5,'Figure 2. Higher-L zoom of the rejected candidate. The negative direction weakens with increasing L but persists.')

# terminology
h1('2. Terminology: what exactly was meant by "closed"')
p('The project accumulated several distinct closure layers. Keeping them separate is essential for later work.')
h2('2.1 Slot closure: 41/41')
p('The common unreduced even-parity coefficient language contains exactly 41 coefficient slots:')
eq('a1...a9  (9)  +  b1...b5  (5)  +  c1...c6  (6)  +  d1...d4  (4)  +  e1...e4  (4)  +  v1...v13  (13)  =  41.')
p('"41/41 closed" means that complete unreduced coefficient tables exist for selected local/sector representatives in this common basis. It does not by itself prove that all radial sectors belong to one identical holonomic action member, nor that the reduced finite-l kinetic matrix is positive everywhere.')
h2('2.2 Constructive same-action closure')
p('This is an existence statement. The outer and inner overlap bands possess enough background-constrained and background-null action-jet freedom to connect the local representatives by smooth C-infinity handovers while retaining the frozen P5 geometry and nonzero common constraint pivots.')
h2('2.3 Principal global closure')
p('A separate center-to-infinity characteristic/principal stream exists. It proves that one can select a globally healthy principal-symbol representation. It is not automatically the same basis, same lower-order member, or same action representative as an independently built finite-l mass/potential block.')
h2('2.4 Finite-l global same-action closure')
p('This is the strongest level required before a coupled QNM spectrum can be called physical: one holonomic center-to-infinity action member, with both handovers directly evaluated from that member, reduced with all profile derivatives retained, and passing K>0 plus radial and angular gradient gates at finite l. This final level is not yet passed.')

# theory architecture
h1('3. Frozen P5 geometry and theory architecture')
h2('3.1 Background variables')
p('The P5 geometry is source-locked in the segment variable Xi. With r_s=1, C=r_s/r=1/x and phi=Xi(C). The susceptibility chi=dXi/dC obeys the autonomous P5 closure dXi/dC=(1-Xi)^(3/2) P5(Xi). The temporal metric is f=(1+Xi)^(-2); h is fixed by the P5 radial closure.')
bullet('Outer unstable light ring: C=2/3, x=1.5.')
bullet('Inner stable light ring: C approximately 0.706135, x approximately 1.41616064.')
bullet('Regular horizonless center.')
p('The historical Hermite middle bridge is provenance only and is not the P5 field equation.')
h2('3.2 Exact kappa identity and the 00/11 basis issue')
eq('phi_r = -C^2 chi,     kappa = h phi_r^2 = h C^4 chi^2 = -2 X.')
p('For the metric background block (f2,f2X), the E00/E11 Jacobian has the triangular form J=[[-A,0],[-A,-B]] with A=r^2 f and B=A kappa. Merely swapping E00 and E11 cannot improve singular values. The algebraically natural basis is instead E0=E00 and EDelta=E11-E00, in which the Jacobian is diagonal: diag(-A,-A kappa).')
p('This is the precise version of the earlier 00/11 rotation/inversion intuition: the useful operation is not a blind row swap but a regime-aware basis change. The same kappa=-2X also controls the rank of the simple SVT Hessian control map. Near the exact center kappa->0, so the scalar coordinate chart must be replaced by the separate regular Taylor-center branch rather than divided through by kappa.')
h2('3.3 Hybrid covariant architecture')
eq('C_HSVT = C_MH + ( C_SVT - C_shared ).')
p('This addition must be performed on the common unreduced quadratic action before eliminating H0,H1,H2,h1,dA0,dA1 or the auxiliary vector variable V. Constraint elimination is nonlinear in the coefficient set, so Reduce(A+B) is not equal to Reduce(A)+Reduce(B) in general.')
p('The physical architecture is therefore Horndeski exterior -> genuine U(1)-SVT light-ring lobe -> Horndeski core, with the genuine-SVT lobe containing both light rings.')

# regions and principal export
h1('4. Established local principal witnesses and the global canonical principal stream')
h2('4.1 Central genuine-SVT reference')
p('The central exact Zhang-Kase reference on 0.61 <= u <= 0.715 is the strongest local even-sector witness. Regressed characteristic minima include:')
table(['Quantity','Certified minimum / range'],[
('K_even','about 10.5415 > 0'),('c_r,T^2','>= 1.00507'),('c_r,S^2','>= 0.903745'),('c_r,V^2','>= 9.45887'),('c_Omega,-^2','>= 84.912'),('c_Omega,+^2','>= 643.142')],widths=[6,8],font=8)
h2('4.2 Global principal export')
g=pd.read_csv(B/'ssz_p5_F2A_GLOBAL_CANONICAL_KRG_PRINCIPAL_CINF_FINAL_2026-09-15(2).csv')
reg=g.groupby('region').size()
rows=[(idx,int(val)) for idx,val in reg.items()]
table(['Region','Rows'],rows,widths=[11,3],font=8)
mins=g[['f','h','cT2','cS2','cV2']].min()
p(f"The 45,617-point principal export has min f={mins.f:.6g}, min h={mins.h:.9g}, min cT^2={mins.cT2:.9g}, min cS^2={mins.cS2:.9g}, min cV^2={mins.cV2:.9g}. All stored R_ij are exactly zero in this canonical export.")
callout('Important limitation', 'This global K/R/G file is a principal-characteristic witness. It must not be silently spliced to an independently constructed finite-l M block and relabeled as a single same-action physical KRGM operator.', 'FCE4D6')

# same action handovers
h1('5. Same-action handovers and background re-solving')
h2('5.1 Outer handover')
p('The explicit re-solved outer background band covers approximately u=0.551523087 to 0.61. The flat partition is genuinely flat at the Horndeski endpoint: S_SVT=0 through the initial plateau, rises across the overlap, and reaches S_SVT=1 at u=0.61. Typical values are approximately S=0.001 near u=0.575, S=0.065 near 0.58, S=0.501 near 0.59, S=0.935 near 0.60, and 1 at 0.61.')
add_picture(B/'fig_outer_partition.png',15.5,'Figure 3. Stored outer C-infinity partition weights T_H and S_SVT.')
p('The resolved metric/vector residuals are at machine precision: max |E00_resolved| about 1.8e-16, max |E11_resolved| about 3.2e-16, max |J_A| about 1.6e-15.')
h2('5.2 Inner handover')
p('The explicit inner band covers u=0.71...0.715. The actual re-solved residuals are machine zero (E00 about 8.9e-16, E11 about 4.7e-15, J_A about 3.1e-15). A previous diagnostic column E_H+E_SVTg is not the complete residual in this band because partition derivatives and re-solved lower Horndeski jets contribute. Treating that naive sum as a failure was superseded.')
p('The stored genuine-SVT lower values themselves visibly follow the action partition: f2=S_SVT*f2_full and f2X=S_SVT*f2X_full to numerical precision, with f2F approximately S_SVT. This is why a split block must not be reinterpreted as an unsplit full-SVT background action.')

# controls
h1('6. Background-null controls and why local existence survived the handovers')
p('The static background fixes only action data tangent to the one-dimensional background trajectory. Transverse jets remain available for perturbative control. For the genuine-SVT Hessian controls H=(HXX,HXF,HXY,HFF,HFY,HYY), the three radial channels are:')
code('Delta v1 = HFF - 4 kappa HFY + 4 kappa^2 HYY\nDelta v4 = -HXF + 2 kappa HXY - 4 Fbg HFY + 8 Fbg kappa HYY\nDelta c2 = (kappa/2) HXX - Fbg HXF + 6 Fbg kappa HXY - 4 Fbg^2 HFY + 16 Fbg^2 kappa HYY')
p('In the selector gauge HXY=HFY=HYY=0, the exact right inverse is:')
eq('HFF = Delta v1,    HXF = -Delta v4,    HXX = (2/kappa)(Delta c2 - Fbg Delta v4).')
eq('det M_rad = kappa/2.')
p('Certified kappa minima are approximately 0.098838 in the outer overlap, 0.207862 in the central lobe, and 0.632753 in the inner overlap, so this control map does not lose rank there. Independent unused Hessian directions remain available for angular-channel regulation. The Horndeski principal-control map with targets (F,H,a1,c4,c2) was likewise found full rank through the onset-to-core interval.')

# derivatives
h1('7. Radial derivative service, holonomicity, and branch corrections')
h2('7.1 JET9D8')
p('The accepted radial derivative service is a local polynomial jet on the nonuniform r=x grid with window=9 and polynomial degree=8. This replaced spline-based derivative experiments for production finite-l work.')
bullet('v12=-v6/(2h) regression: absolute residual around 1e-14 or better.')
bullet('v13: median scaled relative error around 6e-6 on the central regression; edge maxima are derivative/rounding sensitive.')
bullet('e4: isolated median at the few-1e-7 scale; end-to-end regression is also strong.')
p('The previously selected plus branch v12=+v6/(2h) is superseded. Pure Horndeski has v6=0, so the sign was invisible there, but it matters in genuine-SVT and hybrid bands.')
h2('7.2 Correct a5')
eq("a5 = a2' - a1'' - (A0' v4 / 2)' + A0' v5 / 2.")
p("For pure Horndeski A0'=0, this reduces to a5=a2'-a1''. A previous builder omitted -a1'' and is not authoritative.")
h2('7.3 Holonomic on-curve partials')
p("Two Appendix expressions exposed the same practical point: in the selected one-dimensional action representative, partial derivatives appearing in d3 are represented by the holonomic on-curve derivatives, e.g. partial_phi v6=(dv6/dr)/phi_r and partial_phi a4=a4'/phi_r. Using naive explicit partials generated false mismatches. With the holonomic convention, d3 regression closes to numerical/JET accuracy.")
h2('7.4 Central action holonomicity')
p('Direct chain-rule checks such as df3/dr=f3_phi phi_r+f3_X X_r and the corresponding relations for f3_X, f4, f4_X, f4_XX, f2 and f2_X close in the interior, usually around 1e-6 or better and much tighter for several upper jets. Therefore the temporary hypothesis that the central action data were fundamentally nonholonomic was rejected.')

# coefficient emitters
h1('8. Unreduced coefficient emitters')
h2('8.1 Complete Zhang-Kase Appendix-A 41-slot emitter')
p('A complete ZK 41-slot emitter was reconstructed from action jets and frozen P5 geometry. It evaluates a1...a9, b1...b5, c1...c6, d1...d4, e1...e4 and v1...v13 using JET9D8 for radial derivative terms and the selected lower-order member when supplied.')
zk=pd.read_csv(B/'ssz_p5_ZK_APPENDIX_A_EMITTER_OUTER_REGRESSION_2026-09-16.csv')
sel=['a2','c2','d3','e4','v1','v4','v6','v12','v13']
rows=[]
for s in sel:
    rr=zk[zk.slot==s].iloc[0]
    rows.append((s,f'{rr.median_scaled_rel:.3g}',f'{rr.p95_scaled_rel:.3g}',f'{rr.max_scaled_rel:.3g}'))
table(['Slot','median scaled rel','p95','max (edge sensitive)'],rows,font=7.8)
p('The relatively large edge maxima for d3, e4 and v13 are consistent with known high-derivative edge sensitivity and archived rounded inputs; the interior median/p95 values are the meaningful regression measure.')
h2('8.2 Variable-G4(phi), G5=0 Maxwell-Horndeski emitter')
p('The old strong-H principal action dataset was recovered with direct columns for G4(phi), G4_phi, G3_X, F=G=H, a1,a4,c4,mu and P1 over u approximately 0.551523...0.71. This made it possible to reconstruct the correct luminal variable-G4 Horndeski member rather than inventing it backward from a coefficient table.')
mh=pd.read_csv(B/'ssz_p5_MH_G4PHI_EMITTER_STRONGH_REGRESSION_2026-09-16.csv')
sel=['a1','a2','c2','c4','c6','d3','e1','e4','mu','K_scalar']
rows=[]
for s in sel:
    rr=mh[mh.slot==s].iloc[0]
    rows.append((s,f'{rr["median"]:.3g}',f'{rr.p95:.3g}',f'{rr["max"]:.3g}'))
table(['Slot/quantity','median','p95','max'],rows,font=7.8)
p('The Sep-16 parser fixes to c6, d3 and e4 reduce those expressions to numerical/JET-level regression. The emitter is therefore no longer an uncertainty in the final failure diagnosis.')

# constraints
h1('9. Complete common even-parity constraint elimination')
p('Before reduction the common gauge-fixed fields are (H0,H1,H2,h1,dphi,dA0,dA1,V). The physical basis is chi=(psi,dphi,V)^T with')
eq('psi = H2 + (L a4/a3) h1 + (a1/a3) dphi\' .')
p('Define p=a1/a3, q=L a4/a3 and')
code('B    = a2 - v2 v4/(2 v1)\nCphi = a5 + L a6 - v2 v5/(2 v1)\nCH2  = a7 + L a8 - v2 v3/(2 v1)\nCh1  = L (a9 - v2 v6/(2 v1))\nBeff = B - a3 p\' - CH2 p\nDh1  = Ch1 - a3 q\' - CH2 q')
eq("h1 = -(a3 psi' + Beff dphi' + Cphi dphi + CH2 psi + v2 V) / Dh1.")
p("The remaining two algebraic time-derivative constraints use DeltaV=4 b1 v10-v11^2. Finally dA0=(v1 V)'/(L v9) -(v8 h1+v12 H2+v13 dphi)/(2v9). The full JET9D8 implementation exports linear profile maps for h1,H2,H1,dA1,dA0 in y=(psi,dphi,V) and y'.")
p('On the central 4000-point reference, L=6 gives min |Dh1|=11.987940, min |DeltaV|=3.042884 and min |2 L v9|=0.119951. The dphi\'\' and h1\' higher-derivative cancellations close to zero / about 1e-16.')

# profile reducer and S
h1('10. Profile-aware reducer and the missing S block')
p('The decisive structural correction on Sep-16 was to retain the antisymmetric first-order radial block. Zhang-Kase Eq. 4.25 has the reduced form')
eq("L_even^(2) = dot(Y)^T K dot(Y) + Y'^T G Y' + Y'^T S Y + Y^T M Y,    S^T=-S.")
p('In the project convention used by the profile reducer, after applying the formal differential transformation T and integrating by parts:')
code("K = -P20/2\nG = +P02/2        (project action uses -chi'^T G chi')\nR = -P11/2 = 0\nS = 1/2 (P02' - P01)\nM = -1/2 (P00 + S')")
p('Earlier attempts to force the antisymmetric first-order radial contribution into M caused an apparent lower-order mismatch. Once S was retained explicitly, the reducer high-L mass regression passed.')
h2('10.1 Structural reducer audit')
p('On the central reference there are no differential operator terms above second radial order, the H0 quadratic residual is about 4.55e-13, and the internal principal K comparison is typically at 1e-12 scale or better. The project has R=0 in the final profile convention.')

# mass formulas
h1('11. High-L mass/potential closure')
p('The direct-action high-L mass formulas serve as an independent regression target for the profile reducer. The corrected leading terms are:')
code('M11 = -r^2 v6^2/(4 v1)\nM13 = -r v6/2\nM33 = -v1\nM22 = e4 + 2 h c4 a6/a4 - m5_minus^2/(4 a4^2 v9)')
code('z1       = (a4 + r a9 - 0.5 r A0prime v6) v6\nm1_plus  = +r a4 v8 + z1          # c4 sector\nm1_minus = -r a4 v8 + z1          # dA0/v9 Schur sector\nm2       = 2 c5 v1 + (A0prime v1 + 0.5 phiprime v4) v6\nm3_plus  = a4 v1\' + 0.5 A0prime v1 v6\nm4       = 2 A0prime v1 + phiprime v4\nm5_minus = a4 v13 - a6 v6')
code("M12 = -r d3/2 - h c4 z1/(a4 v6) - r h a6 m2/(2 a4 v1) + r v5 v6/(4 v1)\n      - m1_minus m5_minus/(4 a4^2 v9)\n      + (1/4) d/dr [2 r h c4 + r(2 d2 v1-v4 v6)/(2 v1) + r v6 m5_minus/(2 a4 v9)]\n\nM23 = h A0prime c4 v1/a4 + v5/2 + m3_plus m5_minus/(2 a4^2 v9)\n      - h a6 m4/(2 a4) - (1/4) d/dr [v4 + v1 m5_minus/(a4 v9)]")
p('The sign in m5_minus and the distinction m1_plus versus m1_minus are essential. Spline-based mass/angular scans that violated these identities were superseded.')

# 41 stream diagnostics / representative mismatch
h1('12. Why the old 41/41 tables could not simply be concatenated')
p('A diagnostic global 41-stream concatenation was intentionally attempted and rejected. It exposed two separate issues:')
numbered('The archived outer 37/41 table was a repaired pure-SVT branch, not the complete same-action Horndeski+SVT handover. Its coefficients became nonzero before the stored S_SVT partition actually turned on.')
numbered('At u=0.61, the re-solved outer repaired-SVT endpoint and the Sep-13 central exact SVT table are different action-jet representatives on the same background. Their upper/background values f3,f4,N4 are close, but f3X, f2, f2X and Hessian jets are drastically different.')
p('This is not a contradiction: background-constrained jets and background-null transverse jets are not unique. It does mean that a direct coefficient concatenation is not a proof of one same action.')
rows=[]
jumps=pd.read_csv(B/'ssz_p5_OUTER_FINAL_vs_CENTRAL_ENDPOINT_SLOT_JUMPS_2026-09-16.csv')
for _,r in jumps.head(10).iterrows(): rows.append((r.slot,f'{r.outer_u061:.6g}',f'{r.central_u061:.6g}',f'{r.scaled_jump:.3g}'))
table(['Slot','outer at u=.61','central at u=.61','scaled jump'],rows,font=7.6)
callout('Interpretation', 'The existing central_outer_matching_buffer in the global principal export is a characteristic/principal-basis matching section. It does not by itself certify a finite-l action-level 41-slot handover between two different transverse representatives.', 'FFF2CC')

# outer final work timeline
h1('13. Sep-16 outer finite-l endgame - what was tried and what each attempt taught us')
h2('13.1 First simple cubic/Einstein-like MH member')
p('A deliberately simple G4=1/2, G5=0 Maxwell-Horndeski member was combined with the genuine-SVT Delta-C block. Its structural reduction was finite, but K became negative. Large scans of the genuine-SVT Hessian controls, including Delta v1 variations far beyond the initial value, barely moved the bad kinetic direction. The problematic eigenvector was nearly scalar. Therefore the failure was not fixable by simply turning the radial f2 Hessian knobs.')
h2('13.2 Control experiment with the archived strong-H carrier')
p('The same genuine-SVT Delta-C contribution combined at coefficient level with the already healthy strong-H principal carrier gave a high-L K with no negative rows in the diagnostic L=1000 run (minimum eigenvalue approximately 1.18e-6). This isolated the missing freedom to the Horndeski principal member, not the ZK emitter or the reducer.')
h2('13.3 Recovery of the actual strong-H action data')
p('The archived file ssz_p5_horndeski_carrier_through_light_rings_to_core_2026-09-12.csv was recovered. It directly stores the variable G4(phi), G4_phi and G3_X principal carrier. The strong-H table has F=G=H and a positive target scalar kinetic margin about 9.7423e-3. This was the correct input for a variable-G4(phi), G5=0 MH emitter.')
h2('13.4 Auxiliary-V recanonicalization after hybrid assembly')
p('A major algebraic trap was discovered: adding two coefficient sets that had each already been auxiliary-V completed does not preserve the H0^2 identity. After summing the common unreduced blocks, the shared auxiliary V sector must be canonicalized once again. The required identity is v7=v2^2/(4v1) in the corresponding canonical route. A pre-recanonicalization residual around 0.684 was therefore an assembly error, not physics. After common recanonicalization, the residual is exactly zero in the final audit.')
h2('13.5 Final variable-G4 outer candidate - important on-shell caveat')
p('The final outer candidate used:')
bullet('the variable-G4(phi), G5=0 healthy MH principal carrier reconstructed from the archived strong-H action data;')
bullet("the re-solved outer geometry and A0' profile;")
bullet('the genuine-SVT Delta-C block assembled before constraint elimination;')
bullet('common auxiliary-V recanonicalization after assembly;')
bullet('v12=-v6/(2h);')
bullet('selected lower-order member v5=c3=e3=0;')
bullet('a5 regenerated holonomically with JET9D8.')
p('All 41 slots are finite on the 2201-point outer grid, and several assembly identities close. However, this construction is not a strict certified same-action/on-shell member: the variable-G4 strong-H principal data were interpolated onto the already re-solved outer background that had been constructed for a different Horndeski representative, and the full background equations were not re-solved after that substitution. The file names containing FINAL_SAME_ACTION therefore overstate the status and should be read as historical run names. The candidate is nevertheless rejected because the symmetric kinetic quadratic form is negative over a large radial set.')

# final fail details
h1('14. Final variable-G4 candidate operator audit - decisive for the candidate, not a theory no-go')
p('For this candidate, the symmetric part of the 3x3 kinetic quadratic form is not positive definite. The audit code explicitly diagonalizes Ks=(K+K^T)/2. This is the correct quadratic-form test even when a small numerical antisymmetric remainder is present. The detailed audit across L is:')
rows=[]
for _,r in audit.iterrows():
    rows.append((int(r.L),f'{r.min_eig_K:.6e}',int(r.negative_K_rows),f'{r.min_radial_c2:.6e}',int(r.negative_c2_rows),f'{r.max_S_sym:.2e}',f'{r.max_R_abs:.1e}'))
table(['L','min eig(K)','negative rows','min radial c^2','negative c^2 rows','max S sym residual','max |R|'],rows,font=7.4)
p('Additional checks: the H0 auxiliary identity and v12 identity close exactly, no forbidden higher-order operator keys appear, R=0, G is symmetric to floating-point accuracy, and S is antisymmetric at about 1e-8. The raw time-principal block P20/K, however, retains a nonzero antisymmetric numerical remainder in this transplanted candidate (for the stored L=6 and L=42 profiles the relative max asymmetry is about 1.6e-4; at L=1000 about 1.2e-6). This does not affect the quadratic form after symmetrization, but it is another reason not to label the candidate a fully clean final member. Constraint pivots remain nonzero.')
h2('14.1 Spatial location of the kinetic failure')
p('Using the stored L=6,42,1000 matrix profiles, the negative region begins progressively farther inward as L increases:')
rows=[('6','u about 0.552613 to 0.610000','minimum near u=0.608139'),('42','u about 0.554474 to 0.610000','minimum near u=0.608139'),('1000','u about 0.566913 to 0.600830','minimum near u=0.583366')]
table(['L','negative-K interval on stored grid','location of most negative point'],rows,font=8)
p('The high-L weakening explains why principal-only or asymptotic diagnostics can look healthy while finite-l remains unacceptable. It also shows why QNM must wait for an explicitly finite-l healthy member.')
h2('14.2 What the failure does and does not mean')
callout('Does mean', 'The Sep-16 variable-G4 transplant candidate must be rejected. Its symmetrized kinetic quadratic form has a negative direction over a large part of the outer band, and it was not independently re-solved on the background after changing the Horndeski principal member.', 'F4CCCC')
callout('Does not mean', 'It does not show that a certified global same-action finite-l member has failed, because this candidate was not fully background re-solved after the variable-G4 substitution. It also does not invalidate the P5 geometry, local strong-H witness, central genuine-SVT light-ring witness, constructive same-action existence theorem, 41/41 coefficient completion, derivative service, constraint maps, ZK emitter, MH emitter, or profile reducer.', 'D9EAD3')

# QNM
h1('15. QNM status and why calculation stops here')
p('A coupled QNM spectrum is not scientifically meaningful at this stage. The last candidate is rejected and no certified center-to-infinity on-shell finite-l member has yet passed the gates. The QNM gate is therefore deliberately closed.')
p('Independent work on the asymptotic exterior remains useful: an outgoing Jost expansion was transported with order-10 error around 2.43e-12. The previously attempted inward shooting of damped modes on the real radial axis was rejected because the growing complementary solution contaminates the integration and candidate roots moved when the outer radius changed.')
p('Acceptable production approaches after a healthy global KRGM member exists include exterior complex scaling, compactified Jost/spectral determinants, or a continued-fraction-type formulation. Any eventual root must be tested for radial resolution, domain/outer-boundary, complex-scaling angle, basis normalization, branch tracking and determinant convergence.')

# superseded
h1('16. Superseded and rejected routes - do not accidentally reuse')
rows=[
('v12=+v6/(2h)','SUPERSEDED','Accepted Appendix branch is v12=-v6/(2h).'),
("a5 without -a1''",'SUPERSEDED','Correct a5 includes -a1 double prime and vector terms.'),
('Spline-based negative mass/angular scans','SUPERSEDED','They violated later action-level high-L identities.'),
('Central member fundamentally nonholonomic','REJECTED','Direct chain-rule tests close with JET9D8.'),
('A2 requires +4h instead of +2h','REJECTED','Came from a manual derivative error.'),
('Naive inner E_H+E_SVTg residual is the full residual','SUPERSEDED','Re-solved handover includes omitted partition and lower-jet terms.'),
('Eq. A2/131 repair member as production member','NOT ACCEPTED','Diagnostic/existence probe only.'),
('Old F3 inner 41/41 coefficient blend','CANDIDATE ONLY','Pending direct action-level Appendix re-evaluation.'),
('Old outer repaired-SVT 37/41 as full handover','REJECTED FOR HANDOVER','It is a repaired pure-SVT branch, not the same-action partition.'),
('Simple G4=1/2 outer candidate','REJECTED','Too restrictive; produced negative scalar kinetic channel.'),
('Variable-G4 transplanted outer candidate','REJECTED CANDIDATE','Negative sym(K); principal member was transplanted without a fresh full background re-solve.'),
('Real-axis inward QNM shooting','REJECTED','Numerically contaminated for damped modes.')]
table(['Route','Status','Reason'],rows,font=7.4)

# future work
h1('17. Correct future continuation - one targeted control problem, not another tuning loop')
p('If the project is resumed, the next task should not be another hand-picked representative. The correct problem is to formulate the remaining Horndeski/SVT transverse jet freedom as a constrained continuation/optimization problem with the finite-l kinetic gate included explicitly.')
p('1. Freeze the already accepted geometry, partition functions, common basis, v12 sign, lower-order convention, JET9D8 service and reducer.')
p('2. Parameterize only genuinely background-null / representative-selection directions. Do not alter jets already fixed by the frozen background equations unless they are re-solved consistently.')
p('3. Use the full-rank Horndeski control family, not the over-restricted G4=constant subclass. Nonredundant principal targets are (F,H,a1,c4,c2).')
p('4. After selecting a new Horndeski principal member, re-solve the complete static background equations for the combined HSVT action before evaluating finite-l stability. Do not transplant principal jets onto a background solved for another member.')
p('5. Enforce exact background residuals and the common constraint pivots simultaneously with a finite-l lower bound K >= epsilon I on selected L gates, beginning with L=6 and L=42.')
p('6. Only after K>0, solve/check radial generalized eigenvalues and the angular large-L symbol. Then regenerate all 41 coefficients directly from the chosen action jets.')
p('7. Build the action-level Outer->Central and Central->Inner matching sections; do not reuse the principal matching buffers as proof of finite-l action continuity.')
p('8. Freeze one global 41/41 stream and rerun the complete profile reducer. Only then open the QNM gate.')
p('The previously suggested 00/11 inversion/rotation idea remains a useful diagnostic only where the background inverse chart changes conditioning. The mathematically justified form is the E00, E11-E00 basis associated with kappa=-2X. It is not a substitute for kinetic stability control.')

# reproducibility
h1('18. Reproducibility cookbook')
h2('18.1 Authoritative scripts as of this checkpoint')
files=[
('ssz_p5_higher_jet_closure_2026-09-16.py','JET9D8 nonuniform radial derivative/jet service.'),
('ssz_p5_holonomic_a5_closure_2026-09-16.py','Accepted a5 holonomic reconstruction.'),
('ssz_hybrid_full_constraint_maps_JET9D8(1).py','Complete algebraic maps for h1,H2,H1,dA1,dA0.'),
('ssz_hybrid_unreduced_even_kernel.py','Common unreduced action architecture / canonical slot basis.'),
('ssz_hybrid_jet_aware_constraint_reducer.py','First-stage generalized-psi reducer.'),
('ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py','Formal profile-aware T^dagger P T reducer producing K,G,S,M.'),
('ssz_p5_F2b_action_derived_highL_mass_closure_2026-09-15.py','Independent high-L mass identities/regression.'),
('ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py','Complete Zhang-Kase Appendix-A 41-slot emitter.'),
('ssz_p5_mh_luminal_g4phi_emitter_JET9D8_2026-09-16.py','Variable-G4(phi), G5=0 Maxwell-Horndeski emitter after parser fixes.'),
('build_outer_resolved_zk_delta_2026-09-16.py','Outer resolved ZK raw and genuine-SVT Delta-C construction.'),
('final_outer_variableG4_2026-09-16.py','Final Sep-16 variable-G4 transplant candidate assembly/audit; historical naming overstates same-action status.')]
table(['File','Role'],files,widths=[9.5,9],font=7.4)
h2('18.2 Critical data products')
files2=[
('ssz_p5_integrable_svt_unreduced_even_coefficients_2026-09-12.csv','Central raw 41-slot ZK coefficient source.'),
('ssz_p5_integrable_full_svt_lobe_ZK_exact_regression_2026-09-13.csv','Central exact principal/stability regression.'),
('ssz_p5_horndeski_carrier_through_light_rings_to_core_2026-09-12.csv','Recovered strong-H action/principal carrier including G4,G4_phi,G3X.'),
('ssz_p5_F2_horndeski_carrier_unreduced_39of41_CORRECTED_2026-09-14.csv','Archived strong-H coefficient regression target.'),
('ssz_p5_F2_outer_same_action_RESOLVED_background_jets_2026-09-15.csv','Machine-zero outer background handover action jets and partition.'),
('ssz_p5_F2_inner_same_action_RESOLVED_v2_background_jets_2026-09-15.csv','Machine-zero inner background handover action jets.'),
('ssz_p5_F2A_GLOBAL_CANONICAL_KRG_PRINCIPAL_CINF_FINAL_2026-09-15(2).csv','45,617-point global principal K/R/G witness.'),
('ssz_p5_OUTER_GENUINE_SVT_DELTA_41of41_2026-09-16.csv','Outer genuine-SVT beyond-shared coefficient contribution.'),
('ssz_p5_MH_G4PHI_EMITTER_STRONGH_REGRESSION_2026-09-16.csv','Final MH emitter regression statistics.'),
('ssz_p5_ZK_APPENDIX_A_EMITTER_OUTER_REGRESSION_2026-09-16.csv','ZK emitter regression statistics.'),
('ssz_p5_OUTER_FINAL_SAME_ACTION_41of41_2026-09-16.csv','Historical filename: variable-G4 transplant candidate, not a freshly background-re-solved same-action member; candidate fails sym(K).'),
('ssz_p5_OUTER_FINAL_SAME_ACTION_OPERATOR_AUDIT_2026-09-16.csv','Finite-l audit of the rejected transplant candidate; decisive only for that candidate.'),
('SSZ_P5_FINAL_STOP_CHECKPOINT_2026-09-16.md','Short frozen stop-state checkpoint.')]
table(['File','Meaning'],files2,widths=[10.5,8],font=7.3)

h2('18.3 Minimal resume sequence')
code('1. Read SSZ_P5_FINAL_STOP_CHECKPOINT_2026-09-16.md\n2. Read this PDF, especially Sections 2, 10, 13-17\n3. Do NOT regenerate 41/41 from scratch\n4. Do NOT revert v12 sign or a5\n5. Do NOT drop S\n6. Do NOT add reduced operators sector-by-sector\n7. Start from final_outer_variableG4_2026-09-16.py and replace only the representative-selection/control layer\n8. Impose K>0 at finite L as a target during the action-jet continuation\n9. Re-solve background-constrained jets\n10. Re-run full profile reducer and all acceptance gates\n11. Only after PASS: construct coupled QNM determinant')

# provenance/bibliography
h1('19. Literature and convention anchors')
p('The two perturbation formalisms used throughout the coefficient work are:')
bullet('C. Zhang and R. Kase, Physical Review D 110, 044047 (2024), arXiv:2404.11910 - even-parity perturbations in U(1)-gauge-invariant scalar-vector-tensor theories; Appendix-A coefficients and reduced K/G/S/M structure.')
bullet('R. Kase and S. Tsujikawa, Physical Review D 107, 104045 (2023), arXiv:2301.10362 - odd/even black-hole perturbations in Maxwell-Horndeski theories; full five-degree-of-freedom static-spherical formalism.')
p('Project conventions that must stay locked when comparing formulas: r_s=1; u=C=r_s/r; x=r/r_s=1/u; scalar phi=Xi; the accepted vector numbering uses the 13-slot Zhang-Kase schema after canonicalizing the Maxwell-Horndeski v1...v9 vector block; the selected lower-order member used in the final outer run has v5=c3=e3=0.')

# final frozen statement
h1('20. Frozen final statement')
callout('What can safely be said as of 16 September 2026', 'P5 has a regular frozen geometry, healthy local principal witnesses in the exterior, genuine-SVT light-ring lobe and core, constructive same-action handover existence, a complete common 41-slot unreduced coefficient language, a validated derivative service, complete constraint maps, validated Zhang-Kase and Maxwell-Horndeski coefficient emitters, and a profile-aware finite-l reducer including the antisymmetric S term. No certified global finite-l same-action member has yet been completed. The last variable-G4 outer candidate is rejected because sym(K) is negative and, more fundamentally, the changed Horndeski principal member was not followed by a fresh full background re-solve. Therefore coupled QNM remains blocked. Future work should formulate a targeted transverse action-jet continuation/control problem with background equations and finite-l K positivity imposed simultaneously.', 'D9EAD3')
p('This report intentionally records both successes and the final failure so that a later continuation does not accidentally erase hard-won distinctions or redo already closed algebra.')

# appendix: exact audit
h1('Appendix A. Exact final operator-audit rows')
p('A.1 Physical kinetic/radial gate and constraint pivots')
rows=[]
for _,r in audit.iterrows():
    rows.append((int(r.L),f'{r.min_eig_K:.10g}',int(r.negative_K_rows),f'{r.min_radial_c2:.10g}',int(r.negative_c2_rows),f'{r.min_abs_Dh1:.10g}',f'{r.min_abs_DeltaV:.10g}',f'{r.min_abs_2Lv9:.10g}'))
table(['L','min eig K','neg K','min c^2','neg c^2','min |Dh1|','min |DeltaV|','min |2Lv9|'],rows,font=7.0)
p('A.2 Structural residuals of the same final runs')
rows=[]
for _,r in audit.iterrows():
    rows.append((int(r.L),f'{r.max_h0_quadratic:.3g}',f'{r.max_P20_asym:.3g}',f'{r.max_P02_asym:.3g}',f'{r.max_P11_asym:.3g}',f'{r.max_S_sym:.3g}',f'{r.max_R_abs:.3g}',f'{r.max_v12_identity:.3g}'))
table(['L','H0 quadratic','P20 asym','P02 asym','P11 asym','S sym residual','max |R|','v12 residual'],rows,font=7.0)

h1('Appendix B. Final emitter regression highlights')
p('Zhang-Kase outer regression: selected median scaled-relative values.')
rows=[]
for s in ['a2','a7','c2','c5','c6','d3','e1','e4','v1','v4','v9','v13']:
    rr=zk[zk.slot==s].iloc[0]; rows.append((s,f'{rr.median_scaled_rel:.8g}',f'{rr.p95_scaled_rel:.8g}',f'{rr.max_scaled_rel:.8g}',f'{rr.max_abs:.8g}'))
table(['slot','median','p95','max scaled','max abs'],rows,font=7)
p('Maxwell-Horndeski variable-G4 strong-H regression: selected values.')
rows=[]
for s in ['a1','a2','a7','b3','c2','c6','d3','e1','e4','mu','K_scalar']:
    rr=mh[mh.slot==s].iloc[0]; rows.append((s,f'{rr["median"]:.8g}',f'{rr.p95:.8g}',f'{rr["max"]:.8g}'))
table(['slot/quantity','median','p95','max'],rows,font=7)

h1('Appendix C. Compact decision log')
rows=[
('Before Sep-16','41/41 local selected tables existed','Keep: slot closure was real.'),
('Sep-16 early','Profile reducer initially treated lower radial first-order term incorrectly','Fixed by explicit antisymmetric S block.'),
('Sep-16','ZK Appendix-A emitter reconstructed','Keep and reuse.'),
('Sep-16','Outer repaired-SVT endpoint found not identical to central Sep-13 action representative','Do not concatenate as same action.'),
('Sep-16','Simple G4=1/2 MH outer member gave negative K','Reject that restricted member only.'),
('Sep-16','Hessian scans failed to move bad scalar kinetic direction','Do not waste time on f2 Hessian-only tuning.'),
('Sep-16','Strong-H diagnostic restored positive high-L K','Horndeski principal freedom is relevant.'),
('Sep-16','Recovered direct variable-G4 strong-H action data','Use as MH action anchor.'),
('Sep-16','c6,d3,e4 parser/holonomic fixes completed','MH emitter regression closed.'),
('Sep-16','Auxiliary-V block found non-additive after separate completion','Always recanonicalize common V after assembly.'),
('Sep-16 final','Variable-G4 principal carrier transplanted into resolved outer grid and reduced','Rejected candidate: negative sym(K), nonzero raw K asymmetry, and no fresh full background re-solve after transplant.'),
('Stop rule','No more ad-hoc tuning in this run','Resume only as constrained control/optimization problem.')]
table(['Stage','Finding','Rule for future'],rows,font=7.1)

# final metadata
p('Generated from the frozen project artifacts present in the working archive on 16 September 2026. Numerical values in this report are intended as reproducibility anchors; the CSV and Python files listed in Section 18 remain the machine-readable source of truth.', style='Small')

doc.save(OUT)
print(OUT)
