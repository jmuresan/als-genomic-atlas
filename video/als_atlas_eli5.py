"""
ELI-5 explainer for the ALS Genomic Atlas — the 11 data categories.

Render:
    manim -qm video/als_atlas_eli5.py ALSAtlasELI5      # 720p30
    manim -qh video/als_atlas_eli5.py ALSAtlasELI5      # 1080p60

Uses only Text (Pango) and vector shapes, so no LaTeX is required.
"""

import numpy as np
from manim import *

BG = "#0d1b2a"
INK = "#e0e1dd"
DIM = "#415a77"


# --------------------------------------------------------------------------
# Small icon builders. Each returns a VGroup centered near the origin; the
# scene normalizes height before placing it, so absolute sizes here are loose.
# --------------------------------------------------------------------------

def make_pill(color, w=1.4, h=0.55, color2=None):
    pill = RoundedRectangle(
        width=w, height=h, corner_radius=h / 2,
        fill_color=color, fill_opacity=1, stroke_color=WHITE, stroke_width=3,
    )
    divider = Line(UP * (h / 2 - 0.03), DOWN * (h / 2 - 0.03), color=WHITE, stroke_width=3)
    if color2 is not None:
        shade = Rectangle(width=w / 2 - 0.04, height=h - 0.06,
                          fill_color=color2, fill_opacity=1, stroke_width=0)
        shade.move_to(pill).shift(RIGHT * w / 4)
        return VGroup(pill, shade, divider)
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
        Text("G", font_size=40, color=RED, weight=BOLD),
        Text("T", font_size=40, color=WHITE),
    ).arrange(RIGHT, buff=0.18)
    glass = Circle(radius=0.42, color="#ffd166", stroke_width=5).move_to(seq[2])
    handle = Line(glass.get_corner(DR),
                  glass.get_corner(DR) + np.array([0.32, -0.32, 0]),
                  color="#ffd166", stroke_width=6)
    return VGroup(seq, glass, handle)


def icon_regulation(color):
    track = RoundedRectangle(width=1.3, height=0.55, corner_radius=0.275,
                             stroke_color=WHITE, stroke_width=3,
                             fill_color="#1b263b", fill_opacity=1)
    knob = Circle(radius=0.22, color=color, fill_opacity=1, stroke_width=0)
    knob.move_to(track.get_right() + LEFT * 0.3)
    label = Text("ON", font_size=22, color=color).next_to(track, UP, buff=0.12)
    return VGroup(track, knob, label)


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


CATEGORIES = [
    (1,  "Gene & Transcript Map",   ["Find the gene's address", "and read its recipe cards."],   "#4cc9f0", icon_mapping),
    (2,  "Variants & Pathogenicity", ["Spot the spelling mistakes —", "and which ones are harmful."], "#f72585", icon_variants),
    (3,  "Regulation & Epigenomics", ["Find the on/off switches", "and the dimmer knobs."],         "#ffd166", icon_regulation),
    (4,  "Expression & Tissues",    ["See which body parts", "actually use the recipe."],          "#06d6a0", icon_expression),
    (5,  "Pathways & Function",     ["Learn the protein's job", "and which team it's on."],         "#b5179e", icon_pathways),
    (6,  "Network Interactions",    ["See which other proteins", "it holds hands with."],           "#4895ef", icon_interactions),
    (7,  "Drugs & Druggability",    ["Check if a medicine", "can grab onto it."],                   "#ef476f", icon_drug),
    (8,  "3D Structure",            ["See the shape it folds into,", "like origami."],              "#fb8500", icon_structure),
    (9,  "Structural Similarity",   ["Find other proteins", "with the same shape."],                "#8ac926", icon_similarity),
    (10, "Matched-Target Drugs",    ["Borrow medicines made", "for those look-alikes."],           "#ff9e00", icon_matched),
    (11, "Repurposing Candidates",  ["Find similar medicines", "we could reuse."],                  "#a0c4ff", icon_repurpose),
]


class ALSAtlasELI5(Scene):
    def construct(self):
        self.camera.background_color = BG
        self.intro()
        self.run_categories()
        self.outro()

    # ----- intro -------------------------------------------------------
    def intro(self):
        title = Text("The ALS Genomic Atlas", font_size=56, weight=BOLD, color=WHITE)
        sub = Text("explained simply", font_size=30, slant=ITALIC, color="#90a4c4")
        sub.next_to(title, DOWN, buff=0.3)
        self.play(Write(title), run_time=1.2)
        self.play(FadeIn(sub, shift=UP * 0.2), run_time=0.6)
        self.wait(1.0)
        self.play(FadeOut(title), FadeOut(sub), run_time=0.6)

        # a motor neuron
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
        neuron = VGroup(dendrites, soma, axon, terminals).move_to(ORIGIN)
        cap1 = Text("ALS slowly damages motor neurons —", font_size=32, color=INK)
        cap2 = Text("the cells that move your muscles.", font_size=32, color=INK)
        VGroup(cap1, cap2).arrange(DOWN, buff=0.18).to_edge(DOWN, buff=0.9)
        self.play(Create(neuron), run_time=1.4)
        self.play(FadeIn(cap1), FadeIn(cap2), run_time=0.6)
        self.wait(1.6)
        self.play(FadeOut(neuron), FadeOut(cap1), FadeOut(cap2), run_time=0.6)

        # genes + the detective framing
        line1 = Text("Genes hold the body's instructions.", font_size=36, color=WHITE)
        line2 = Text("46 genes are linked to ALS.", font_size=36, color="#90a4c4")
        VGroup(line1, line2).arrange(DOWN, buff=0.3)
        self.play(FadeIn(line1, shift=UP * 0.2))
        self.play(FadeIn(line2, shift=UP * 0.2))
        self.wait(1.4)
        self.play(FadeOut(line1), FadeOut(line2), run_time=0.5)

        glass = Circle(radius=0.6, color="#ffd166", stroke_width=6)
        handle = Line(glass.get_corner(DR), glass.get_corner(DR) + np.array([0.5, -0.5, 0]),
                      color="#ffd166", stroke_width=8)
        detective = VGroup(glass, handle).scale(0.9).shift(UP * 0.5)
        framing1 = Text("For each gene, scientists gather 11 kinds of clues", font_size=32, color=INK)
        framing2 = Text("— like a detective building a case.", font_size=32, color="#90a4c4")
        VGroup(framing1, framing2).arrange(DOWN, buff=0.18).to_edge(DOWN, buff=1.1)
        self.play(Create(detective), run_time=0.9)
        self.play(FadeIn(framing1), FadeIn(framing2), run_time=0.6)
        self.wait(1.6)
        self.play(FadeOut(detective), FadeOut(framing1), FadeOut(framing2), run_time=0.6)

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
        for idx, (num, title, eli5, color, icon_fn) in enumerate(CATEGORIES):
            self.set_progress(idx)

            badge = VGroup(
                Circle(radius=0.5, color=color, fill_opacity=1, stroke_width=0),
                Text(str(num), font_size=40, weight=BOLD, color=BG),
            )
            badge[1].move_to(badge[0])
            title_t = Text(title, font_size=38, weight=BOLD, color=WHITE)
            header = VGroup(badge, title_t).arrange(RIGHT, buff=0.35).to_edge(UP, buff=1.0)

            tag = Text(f"Clue {num} of 11", font_size=24, color=color)
            tag.next_to(header, DOWN, buff=0.3)

            icon = icon_fn(color)
            icon.scale_to_fit_height(2.0).move_to(LEFT * 3.3 + DOWN * 0.3)

            eli5_t = VGroup(*[Text(t, font_size=32, color=INK) for t in eli5])
            eli5_t.arrange(DOWN, aligned_edge=LEFT, buff=0.22)
            eli5_t.next_to(icon, RIGHT, buff=1.0).set_y(icon.get_y())

            self.play(FadeIn(badge, scale=0.6), Write(title_t), FadeIn(tag), run_time=0.7)
            self.play(Create(icon), run_time=0.9)
            self.play(*[FadeIn(t, shift=RIGHT * 0.3) for t in eli5_t], run_time=0.7)
            self.wait(2.2)
            self.play(FadeOut(header), FadeOut(tag), FadeOut(icon), FadeOut(eli5_t), run_time=0.45)

    # ----- outro -------------------------------------------------------
    def outro(self):
        self.play(*[
            self.dots[j].animate.set_fill(CATEGORIES[j][3], opacity=1).set(width=0.26)
            for j in range(len(CATEGORIES))
        ], run_time=0.8)

        line1 = Text("11 clues  →  one full picture of a gene.", font_size=38, color=WHITE)
        line1.move_to(UP * 0.6)
        self.play(Write(line1), run_time=1.0)
        self.wait(1.2)

        line2 = Text("× 46 genes  =  the ALS Genomic Atlas", font_size=40, weight=BOLD, color="#4cc9f0")
        line2.next_to(line1, DOWN, buff=0.5)
        self.play(FadeIn(line2, shift=UP * 0.2), run_time=0.8)
        self.wait(1.8)
        self.play(FadeOut(line1), FadeOut(line2), FadeOut(self.dots), run_time=0.8)
        self.wait(0.4)
