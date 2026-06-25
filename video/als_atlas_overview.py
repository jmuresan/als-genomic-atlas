"""
Overview of the ALS Genomic Atlas — the 11 data categories.

Render:
    manim -qm video/als_atlas_overview.py ALSGenomicAtlas      # 720p30
    manim -qh video/als_atlas_overview.py ALSGenomicAtlas      # 1080p60

Uses Text (Pango) and vector shapes, so no LaTeX is required.
The SANDO logo is composited from video/assets/sando_logo_soft.png.
"""

import numpy as np
from manim import *
import os
import subprocess

def get_audio_duration(filename):
    filepath = os.path.join("video", "audio", f"{filename}.mp3")
    if not os.path.exists(filepath):
        return 3.0
    cmd = f"ffprobe -v error -show_entries format=duration -of default=noprint_wrappers=1:nokey=1 {filepath}"
    try:
        res = subprocess.check_output(cmd.split())
        return float(res)
    except Exception:
        return 3.0


BG = "#0d1b2a"
INK = "#e0e1dd"
DIM = "#415a77"
MUTED = "#90a4c4"
LOGO = "video/assets/sando_logo_soft.png"


# --------------------------------------------------------------------------
# Icons. Each returns a VGroup centered near the origin; the scene normalizes
# height before placing it, so absolute sizes here are loose. The icons are
# schematic, not literal depictions.
# --------------------------------------------------------------------------

def make_pill(color, w=1.4, h=0.55):
    pill = RoundedRectangle(
        width=w, height=h, corner_radius=h / 2,
        fill_color=color, fill_opacity=1, stroke_color=WHITE, stroke_width=3,
    )
    divider = Line(UP * (h / 2 - 0.03), DOWN * (h / 2 - 0.03), color=WHITE, stroke_width=3)
    return VGroup(pill, divider)


def icon_mapping(color):
    bar = RoundedRectangle(width=2.4, height=0.4, corner_radius=0.15,
                           stroke_color=WHITE, stroke_width=3,
                           fill_color="#1b263b", fill_opacity=1)
    segs = VGroup(*[
        Rectangle(width=0.45, height=0.26, fill_color=color, fill_opacity=1, stroke_width=0)
        for _ in range(3)
    ]).arrange(RIGHT, buff=0.25).move_to(bar)
    pin_c = Circle(radius=0.2, color=color, fill_opacity=1, stroke_width=0)
    pin_t = Triangle(fill_color=color, fill_opacity=1, stroke_width=0).scale(0.16).rotate(PI)
    pin = VGroup(pin_c, pin_t).arrange(DOWN, buff=-0.05)
    pin.next_to(bar, UP, buff=0.0).align_to(segs[0], LEFT).shift(RIGHT * 0.1)
    return VGroup(bar, segs, pin)


def icon_variants(color):
    seq = VGroup(
        Text("A", font_size=40, color=WHITE),
        Text("C", font_size=40, color=WHITE),
        Text("G", font_size=40, color=color, weight=BOLD),
        Text("T", font_size=40, color=WHITE),
    ).arrange(RIGHT, buff=0.18)
    glass = Circle(radius=0.42, color=WHITE, stroke_width=5).move_to(seq[2])
    handle = Line(glass.get_corner(DR),
                  glass.get_corner(DR) + np.array([0.32, -0.32, 0]),
                  color=WHITE, stroke_width=6)
    return VGroup(seq, glass, handle)


def icon_regulation(color):
    track = RoundedRectangle(width=1.3, height=0.55, corner_radius=0.275,
                             stroke_color=WHITE, stroke_width=3,
                             fill_color="#1b263b", fill_opacity=1)
    knob = Circle(radius=0.22, color=color, fill_opacity=1, stroke_width=0)
    knob.move_to(track.get_right() + LEFT * 0.3)
    return VGroup(track, knob)


def icon_expression(color):
    heights = [0.6, 1.2, 0.85, 1.5]
    bars = VGroup(*[
        Rectangle(width=0.32, height=h, fill_color=color, fill_opacity=1, stroke_width=0)
        for h in heights
    ]).arrange(RIGHT, buff=0.16, aligned_edge=DOWN)
    base = Line(bars.get_corner(DL) + LEFT * 0.12, bars.get_corner(DR) + RIGHT * 0.12,
                color=WHITE, stroke_width=3)
    return VGroup(bars, base)


def icon_pathways(color):
    pts = [LEFT * 1.1 + DOWN * 0.2, LEFT * 0.3 + UP * 0.5,
           RIGHT * 0.4 + DOWN * 0.4, RIGHT * 1.1 + UP * 0.3]
    edges = VGroup(*[Line(pts[i], pts[i + 1], color=WHITE, stroke_width=3) for i in range(3)])
    nodes = VGroup(*[Dot(p, radius=0.14, color=color) for p in pts])
    return VGroup(edges, nodes)


def icon_interactions(color):
    center = Dot(ORIGIN, radius=0.2, color=color)
    edges, nodes = VGroup(), VGroup()
    for ang in (30, 110, 200, 290):
        p = 1.0 * np.array([np.cos(np.deg2rad(ang)), np.sin(np.deg2rad(ang)), 0])
        edges.add(Line(ORIGIN, p, color=WHITE, stroke_width=3))
        nodes.add(Dot(p, radius=0.14, color=color))
    return VGroup(edges, center, nodes)


def icon_drug(color):
    return make_pill(color).rotate(-PI / 6)


def icon_structure(color):
    s = 0.9
    front = Square(side_length=s, color=color, stroke_width=3)
    back = Square(side_length=s, color=color, stroke_width=3).shift(UP * 0.4 + RIGHT * 0.4)
    fv, bv = front.get_vertices(), back.get_vertices()
    conn = VGroup(*[Line(fv[i], bv[i], color=color, stroke_width=3) for i in range(4)])
    return VGroup(back, conn, front)


def icon_similarity(color):
    p1 = RegularPolygon(n=6, color=color, fill_opacity=0.35, stroke_width=3).scale(0.55)
    p2 = p1.copy()
    grp = VGroup(p1, p2).arrange(RIGHT, buff=0.95)
    approx = Text("≈", font_size=54, color=WHITE).move_to(grp.get_center())
    return VGroup(grp, approx)


def icon_matched(color):
    look = RegularPolygon(n=6, color=color, fill_opacity=0.35, stroke_width=3).scale(0.42)
    arrow = Arrow(LEFT * 0.25, RIGHT * 0.25, color=WHITE, buff=0, stroke_width=4)
    pill = make_pill(color, w=1.0, h=0.45)
    return VGroup(look, arrow, pill).arrange(RIGHT, buff=0.25)


def icon_repurpose(color):
    pill1 = make_pill(color, w=0.9, h=0.42)
    pill2 = make_pill("#778da9", w=0.9, h=0.42)
    pills = VGroup(pill1, pill2).arrange(RIGHT, buff=0.8)
    top = CurvedArrow(pill1.get_top() + UP * 0.06, pill2.get_top() + UP * 0.06,
                      angle=-PI / 2, color=WHITE, stroke_width=3, tip_length=0.16)
    bot = CurvedArrow(pill2.get_bottom() + DOWN * 0.06, pill1.get_bottom() + DOWN * 0.06,
                      angle=-PI / 2, color=WHITE, stroke_width=3, tip_length=0.16)
    return VGroup(pills, top, bot)


# (number, short title, two-line factual description, accent color, icon, sources)
CATEGORIES = [
    (1,  "Gene & Transcript Map",   ["The gene's location in the genome", "and its transcript isoforms."],   "#4cc9f0", icon_mapping,     "Ensembl, NCBI, UniProt"),
    (2,  "Variants & Pathogenicity", ["Known DNA variants and their", "assessed clinical significance."],     "#5e8bd6", icon_variants,    "ClinVar, dbSNP, gnomAD, AlphaGenome"),
    (3,  "Regulation & Epigenomics", ["Promoters, enhancers, and", "transcription-factor binding sites."],     "#ffd166", icon_regulation,  "ENCODE, UCSC, JASPAR, UniBind"),
    (4,  "Expression & Tissues",    ["Where the gene is expressed,", "including nervous-system tissue."],       "#06d6a0", icon_expression,  "GTEx, Human Protein Atlas"),
    (5,  "Pathways & Function",     ["Biological pathways, gene functions,", "and protein domains."],            "#b5179e", icon_pathways,    "Reactome, QuickGO, InterPro"),
    (6,  "Network Interactions",    ["Proteins it physically or", "functionally interacts with."],              "#4895ef", icon_interactions, "STRING"),
    (7,  "Drugs & Druggability",    ["Associated drugs, clinical trials,", "and target–disease evidence."], "#ef476f", icon_drug,        "Open Targets, ChEMBL, ClinicalTrials.gov"),
    (8,  "3D Structure",            ["Experimental and predicted", "3D protein structures."],                   "#fb8500", icon_structure,   "AlphaFold, RCSB PDB"),
    (9,  "Structural Similarity",   ["Other proteins with a", "similar 3D structure."],                         "#8ac926", icon_similarity,  "Foldseek"),
    (10, "Matched-Target Drugs",    ["Drugs and trials for those", "structurally similar proteins."],           "#ff9e00", icon_matched,     "Foldseek hits via Open Targets"),
    (11, "Repurposing Candidates",  ["Chemically similar compounds", "flagged for further study."],             "#a0c4ff", icon_repurpose,   "ChEMBL, Open Targets"),
]


class ALSGenomicAtlas(Scene):
    def construct(self):
        self.camera.background_color = BG
        self.intro()
        self.run_categories()
        self.outro()

    # ----- intro -------------------------------------------------------
    def intro(self):
        title = Text("ALS Genomic Atlas", font_size=58, weight=BOLD, color=WHITE)
        logo = ImageMobject(LOGO).scale_to_fit_height(1.1)
        group = Group(title, logo).arrange(DOWN, buff=0.5)
        self.play(Write(title), run_time=1.2)
        self.play(FadeIn(logo), run_time=0.8)
        dur = get_audio_duration("01_intro_title")
        self.add_sound("video/audio/01_intro_title.mp3")
        self.wait(max(0.5, dur - 2.0))
        self.play(FadeOut(title), FadeOut(logo), run_time=0.7)

        # motor neuron schematic
        soma = Circle(radius=0.55, color="#4cc9f0", fill_opacity=0.25, stroke_width=3)
        dendrites = VGroup(*[
            Line(soma.get_center(),
                 soma.get_center() + 1.1 * np.array([np.cos(a), np.sin(a), 0]),
                 color="#4cc9f0", stroke_width=3)
            for a in np.deg2rad([100, 140, 180, 220, 260])
        ])
        axon = Line(soma.get_right(), soma.get_right() + RIGHT * 4.2, color="#4cc9f0", stroke_width=4)
        terminals = VGroup(*[
            Line(axon.get_end(), axon.get_end() + 0.5 * np.array([np.cos(a), np.sin(a), 0]),
                 color="#4cc9f0", stroke_width=3)
            for a in np.deg2rad([30, 0, -30])
        ])
        neuron = VGroup(dendrites, soma, axon, terminals).move_to(UP * 0.4)
        cap = VGroup(
            Text("Amyotrophic lateral sclerosis is a progressive", font_size=30, color=INK),
            Text("neurodegenerative disease that affects motor neurons.", font_size=30, color=INK),
        ).arrange(DOWN, buff=0.16).to_edge(DOWN, buff=1.0)
        dur = get_audio_duration("02_intro_neuron")
        self.add_sound("video/audio/02_intro_neuron.mp3")
        self.play(Create(neuron), run_time=1.4)
        self.play(FadeIn(cap), run_time=0.6)
        self.wait(max(0.5, dur - 2.0))
        self.play(FadeOut(neuron), FadeOut(cap), run_time=0.6)

        framing = VGroup(
            Text("The atlas compiles molecular data for 46 genes", font_size=34, color=WHITE),
            Text("associated with ALS, organized into 11 categories.", font_size=34, color=WHITE),
        ).arrange(DOWN, buff=0.2)
        dur = get_audio_duration("03_intro_categories")
        self.add_sound("video/audio/03_intro_categories.mp3")
        self.play(FadeIn(framing[0], shift=UP * 0.2))
        self.play(FadeIn(framing[1], shift=UP * 0.2))
        self.wait(max(0.5, dur - 0.4))
        self.play(FadeOut(framing), run_time=0.6)

    # ----- progress dots ----------------------------------------------
    def build_progress(self):
        self.dots = VGroup(*[Dot(radius=0.1, color=DIM) for _ in CATEGORIES])
        self.dots.arrange(RIGHT, buff=0.35).to_edge(DOWN, buff=0.5)
        self.play(FadeIn(self.dots), run_time=0.5)

    def set_progress(self, idx):
        for j, d in enumerate(self.dots):
            if j == idx:
                d.set_fill(CATEGORIES[j][3], opacity=1).set_stroke(width=0)
                d.width = 0.30
            elif j < idx:
                d.set_fill(CATEGORIES[j][3], opacity=0.7).set_stroke(width=0)
                d.width = 0.18
            else:
                d.set_fill(DIM, opacity=1).set_stroke(width=0)
                d.width = 0.18

    # ----- the 11 categories ------------------------------------------
    def run_categories(self):
        self.build_progress()
        for idx, (num, title, desc, color, icon_fn, sources) in enumerate(CATEGORIES):
            self.set_progress(idx)

            badge = VGroup(
                Circle(radius=0.5, color=color, fill_opacity=1, stroke_width=0),
                Text(str(num), font_size=40, weight=BOLD, color=BG),
            )
            badge[1].move_to(badge[0])
            title_t = Text(title, font_size=38, weight=BOLD, color=WHITE)
            header = VGroup(badge, title_t).arrange(RIGHT, buff=0.35).to_edge(UP, buff=1.0)

            tag = Text(f"Category {num} / 11", font_size=24, color=color)
            tag.next_to(header, DOWN, buff=0.3)

            icon = icon_fn(color)
            icon.scale_to_fit_height(2.0).move_to(LEFT * 3.3 + DOWN * 0.2)

            desc_t = VGroup(*[Text(t, font_size=31, color=INK) for t in desc])
            desc_t.arrange(DOWN, aligned_edge=LEFT, buff=0.22)
            desc_t.next_to(icon, RIGHT, buff=1.0).set_y(icon.get_y() + 0.25)

            src_t = Text(f"Source: {sources}", font_size=22, color=MUTED)
            src_t.next_to(desc_t, DOWN, aligned_edge=LEFT, buff=0.45)

            filename = f"{(num+3):02d}_cat_{num}"
            dur = get_audio_duration(filename)
            self.add_sound(f"video/audio/{filename}.mp3")
            self.play(FadeIn(badge, scale=0.6), Write(title_t), FadeIn(tag), run_time=0.7)
            self.play(Create(icon), run_time=0.9)
            self.play(*[FadeIn(t, shift=RIGHT * 0.3) for t in desc_t], run_time=0.6)
            self.play(FadeIn(src_t), run_time=0.4)
            self.wait(max(0.5, dur - 2.6))
            self.play(FadeOut(header), FadeOut(tag), FadeOut(icon),
                      FadeOut(desc_t), FadeOut(src_t), run_time=0.45)

    # ----- outro -------------------------------------------------------
    def outro(self):
        self.play(*[
            self.dots[j].animate.set_fill(CATEGORIES[j][3], opacity=1).set(width=0.24)
            for j in range(len(CATEGORIES))
        ], run_time=0.7)

        closing = VGroup(
            Text("These 11 categories are compiled for each", font_size=34, color=WHITE),
            Text("of the 46 ALS-associated genes in the atlas.", font_size=34, color=WHITE),
        ).arrange(DOWN, buff=0.2)
        dur = get_audio_duration("15_outro_summary")
        self.add_sound("video/audio/15_outro_summary.mp3")
        self.play(FadeIn(closing), run_time=0.8)
        self.wait(max(0.5, dur - 0.8))
        self.play(FadeOut(closing), FadeOut(self.dots), run_time=0.6)

        credit = VGroup(
            Text("All data is drawn from public biological databases,", font_size=26, color=MUTED),
            Text("including Ensembl, UniProt, ClinVar, gnomAD, GTEx, STRING,", font_size=26, color=MUTED),
            Text("Reactome, Open Targets, AlphaFold, the RCSB PDB,", font_size=26, color=MUTED),
            Text("Foldseek, and ChEMBL.", font_size=26, color=MUTED),
        ).arrange(DOWN, buff=0.18)
        dur = get_audio_duration("16_outro_credits")
        self.add_sound("video/audio/16_outro_credits.mp3")
        self.play(FadeIn(credit), run_time=0.8)
        self.wait(max(0.5, dur - 0.8))
        self.play(FadeOut(credit), run_time=0.6)

        logo = ImageMobject(LOGO).scale_to_fit_height(2.4)
        self.play(FadeIn(logo), run_time=0.9)
        self.wait(1.8)
        self.play(FadeOut(logo), run_time=0.8)
        self.wait(0.3)
