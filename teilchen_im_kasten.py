# -*- coding: utf-8 -*-
"""
Vom Dünger zur Quantenwelt – Arbeitskreis
Block: Quantenmechanik – Teilchen im Kasten

Diese App ersetzt den Impulsvortrag. Links durchklicken, rechts erklären
und mitspielen lassen.

Start:  python -m streamlit run teilchen_im_kasten.py
Braucht: streamlit, numpy, matplotlib
"""

import numpy as np
import matplotlib.pyplot as plt
import streamlit as st

# ==================================================================
# Konstanten
# ==================================================================
H = 6.62607015e-34       # Planck-Konstante [J*s]
HBAR = H / (2 * np.pi)
C = 299792458.0          # Lichtgeschwindigkeit [m/s]
M_E = 9.1093837015e-31   # Elektronenmasse [kg]
EV = 1.602176634e-19     # 1 Elektronenvolt in Joule

BLAU = "#1b6ca8"
ORANGE = "#d95f02"
GRAU = "#c9c9c9"

st.set_page_config(page_title="Teilchen im Kasten", layout="centered")


# ==================================================================
# Hilfsfunktionen
# ==================================================================
def energie(n, L, m=M_E):
    """Energie des n-ten Zustands im Kasten der Breite L."""
    return (n ** 2 * H ** 2) / (8 * m * L ** 2)


def psi(n, L, x):
    """Wellenfunktion des n-ten Zustands."""
    return np.sqrt(2 / L) * np.sin(n * np.pi * x / L)


def wellenlaenge_zu_farbe(wl_nm):
    """Grobe Umrechnung sichtbare Wellenlänge -> RGB-Hex. Nur zur Anschauung."""
    wl = wl_nm
    if wl < 380 or wl > 750:
        return "#444444"
    if wl < 440:
        r, g, b = -(wl - 440) / 60, 0.0, 1.0
    elif wl < 490:
        r, g, b = 0.0, (wl - 440) / 50, 1.0
    elif wl < 510:
        r, g, b = 0.0, 1.0, -(wl - 510) / 20
    elif wl < 580:
        r, g, b = (wl - 510) / 70, 1.0, 0.0
    elif wl < 645:
        r, g, b = 1.0, -(wl - 645) / 65, 0.0
    else:
        r, g, b = 1.0, 0.0, 0.0
    return "#%02x%02x%02x" % (int(r * 255), int(g * 255), int(b * 255))


def plot_zustand(n, L_nm, modus="welle"):
    L = L_nm * 1e-9
    x = np.linspace(0, L, 600)
    y = psi(n, L, x)
    x_nm = x * 1e9

    fig, ax = plt.subplots(figsize=(7, 3.2))
    if modus == "wahrscheinlichkeit":
        y = y ** 2
        ax.fill_between(x_nm, y, alpha=0.4, color=ORANGE)
        ax.plot(x_nm, y, color=ORANGE, lw=2)
        ax.set_ylabel("Wahrscheinlichkeit")
        ax.set_ylim(bottom=0)
    else:
        ax.plot(x_nm, y, color=BLAU, lw=2.5)
        ax.axhline(0, color="grey", lw=0.8)
        ax.set_ylabel("Welle")
    ax.axvline(0, color="black", lw=5)
    ax.axvline(L_nm, color="black", lw=5)
    ax.set_xlabel("Ort im Kasten [Nanometer]")
    ax.set_xlim(-0.05 * L_nm, 1.05 * L_nm)
    ax.set_yticks([])
    ax.spines[["top", "right", "left"]].set_visible(False)
    return fig


# ==================================================================
# Navigation
# ==================================================================
KAPITEL = [
    "0 · Worum geht es hier?",
    "1 · Die Gitarrensaite",
    "2 · Das Teilchen im Kasten",
    "3 · Die Energieleiter",
    "4 · Farbe aus dem Sprung",
    "5 · Wo ist das Teilchen?",
    "6 · Warum merkst du nichts?",
    "7 · Diskussion",
]

st.sidebar.title("Vom Dünger zur Quantenwelt")
st.sidebar.caption("Block Quantenmechanik")
kapitel = st.sidebar.radio("Kapitel", KAPITEL, label_visibility="collapsed")
st.sidebar.divider()
st.sidebar.caption(
    "Diese App braucht kein Vorwissen. "
    "Alle Fachwörter stehen in Kapitel 7."
)


# ==================================================================
# Kapitel 0
# ==================================================================
if kapitel == KAPITEL[0]:
    st.title("Worum geht es hier?")
    st.markdown(
        """
Um 1900 war die Physik ziemlich zufrieden mit sich. Man konnte Planeten
berechnen, Dampfmaschinen bauen, Brücken auslegen. Die Welt schien im Prinzip
verstanden – es fehlten nur noch ein paar Nachkommastellen.

Dann kamen ein paar Experimente, die einfach **nicht passten**. Und zwar nicht
ein bisschen daneben, sondern grundsätzlich.

Das bekannteste Problem: Ein Atom besteht aus einem Kern und Elektronen, die
um ihn herum sind. Nach der klassischen Physik müsste ein kreisendes Elektron
ständig Energie abstrahlen, langsamer werden und **innerhalb von Sekundenbruchteilen in den Kern stürzen**.

Tut es aber nicht. Du bist der Beweis: Du bestehst aus Atomen, und die sind
seit Milliarden Jahren stabil.

Die klassische Physik konnte nicht erklären, warum es dich gibt.
"""
    )
    st.info(
        "**Die Antwort**, die daraus entstand, heißt Quantenmechanik. "
        "Sie ist berüchtigt dafür, unverständlich zu sein. "
        "Ist sie aber nicht, zumindest nicht der Kern davon. "
        "Den bauen wir in den nächsten Minuten auf."
    )
    st.markdown(
        """
Unser Werkzeug dafür ist das einfachste Beispiel, das es gibt: **ein Teilchen,
das in einem Kasten eingesperrt ist.** Mehr nicht. Daran kann man fast alles
zeigen, was an der Quantenmechanik neu und verstörend war.
"""
    )
    st.caption("→ Weiter zu Kapitel 1 in der Navigation links.")


# ==================================================================
# Kapitel 1
# ==================================================================
elif kapitel == KAPITEL[1]:
    st.title("Die Gitarrensaite")
    st.markdown(
        """
Bevor wir zu Teilchen kommen, etwas Vertrautes.

Eine Gitarrensaite ist an **beiden Enden festgeklemmt**. Genau das ist der
entscheidende Punkt. Wenn du sie zupfst, kann sie nicht irgendwie schwingen,
sie kann nur so schwingen, dass sie an den Enden stillsteht.

Es passen also nur bestimmte Wellen hinein:
"""
    )

    saite = st.slider("Wie viele Bäuche soll die Saite haben?", 1, 6, 1)

    x = np.linspace(0, 1, 500)
    fig, ax = plt.subplots(figsize=(7, 2.6))
    for phase in np.linspace(-1, 1, 7):
        ax.plot(x, phase * np.sin(saite * np.pi * x), color=BLAU, alpha=0.25, lw=1)
    ax.plot(x, np.sin(saite * np.pi * x), color=BLAU, lw=2.5)
    ax.plot(x, -np.sin(saite * np.pi * x), color=BLAU, lw=2.5)
    ax.axvline(0, color="black", lw=6)
    ax.axvline(1, color="black", lw=6)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines[:].set_visible(False)
    st.pyplot(fig)

    st.markdown(
        f"""
Bei **{saite} {'Bauch' if saite == 1 else 'Bäuchen'}** klingt die Saite in einem
ganz bestimmten Ton. Dazwischen gibt es nichts. Es gibt keinen Ton „zwischen"
dem ersten und dem zweiten Oberton – jedenfalls nicht auf dieser Saite mit
dieser Länge.

Und jetzt die zwei Beobachtungen, auf die alles Weitere aufbaut:
"""
    )
    st.success(
        "**1. Eingesperrt sein erzeugt Auswahl.** Während eine freie Welle im Raum jede beliebige Form annehmen kann, kann eine eingeklemmte Welle das nicht mehr."
    )
    st.success(
        "**2. Die Auswahl ist abzählbar.** Du kannst die erlaubten Schwingungen "
        "durchnummerieren: die erste, die zweite, die dritte. Keine "
        "Zwischenstufen."
    )
    st.markdown(
        """
Das ist keine geheimnisvolle Physik. Das ist Handwerk – jeder Instrumentenbauer
weiß das seit Jahrhunderten.

**Die einzige Behauptung der Quantenmechanik lautet: Materie macht das auch.**
"""
    )


# ==================================================================
# Kapitel 2
# ==================================================================
elif kapitel == KAPITEL[2]:
    st.title("Das Teilchen im Kasten")
    st.markdown(
        """
Jetzt sperren wir statt einer Saite ein **Elektron** ein. In einen winzigen
Kasten, aus dem es nicht heraus kann.

Der überraschende Befund der Quantenmechanik: Ein eingesperrtes Teilchen
verhält sich wie eine eingeklemmte Saite. Es hat eine **Welle**, und diese Welle
muss an den Wänden auf null gehen.
"""
    )

    c1, c2 = st.columns(2)
    with c1:
        n = st.slider("Quantenzahl n", 1, 8, 1)
    with c2:
        L_nm = st.slider("Kastenbreite L [Nanometer]", 0.2, 3.0, 1.0, 0.1)

    st.pyplot(plot_zustand(n, L_nm, "welle"))

    L = L_nm * 1e-9
    E_n = energie(n, L)

    k1, k2, k3 = st.columns(3)
    k1.metric("Energie", f"{E_n / EV:.2f} eV")
    k2.metric("Bäuche", f"{n}")
    k3.metric("Nullstellen innen", f"{n - 1}")

    st.divider()
    st.subheader("Was ist die Quantenzahl n?")
    st.markdown(
        """
**n ist einfach eine Hausnummer.** Sie zählt durch, um welche der erlaubten
Wellen es sich handelt.

- n = 1 → die einfachste mögliche Welle, ein Bauch. Der ruhigste Zustand.
- n = 2 → zwei Bäuche, mehr Zappeln, mehr Energie.
- n = 3 → drei Bäuche. Und so weiter.

n ist immer eine **ganze Zahl**. Es gibt kein n = 1,5.
"""
    )
    st.info(
        "**Merksatz für die Runde:** Eine Quantenzahl ist keine Messgröße, "
        "sondern eine Abzählung. Sie sagt nicht *wie viel*, sondern *welcher*."
    )

    st.subheader("Und was ist diese „Welle\" überhaupt?")
    st.markdown(
        """
**Es ist keine Welle aus irgendeinem Stoff.** Das Elektron wackelt nicht auf und
ab wie ein Seil.

Diese Welle heißt **Wellenfunktion**, geschrieben ψ (griechisch „psi").
Sie ist eine Rechengröße. Ihre Bedeutung ist:

> Wo die Welle groß ist, findest du das Teilchen wahrscheinlich.
> Wo sie null ist, findest du es nie.

Das ist die eigentliche Zumutung der Theorie, und darum geht es in Kapitel 5.
"""
    )
    st.markdown(
        "**Probier es aus:** Zieh die Kastenbreite kleiner. Was passiert mit der "
        "Energie? Warum wohl? (Antwort in Kapitel 3.)"
    )


# ==================================================================
# Kapitel 3
# ==================================================================
elif kapitel == KAPITEL[3]:
    st.title("Die Energieleiter")
    st.markdown(
        """
Jede erlaubte Welle gehört zu einer bestimmten Energie. Da es nur bestimmte
Wellen gibt, gibt es auch **nur bestimmte Energien**.
"""
    )

    L_nm = st.slider("Kastenbreite L [Nanometer]", 0.2, 3.0, 1.0, 0.1, key="l3")
    n_sel = st.slider("Markierter Zustand n", 1, 8, 1, key="n3")
    L = L_nm * 1e-9

    fig, ax = plt.subplots(figsize=(7, 3.6))
    for k in range(1, 9):
        E_k = energie(k, L) / EV
        farbe = ORANGE if k == n_sel else GRAU
        ax.hlines(E_k, 0, 1, color=farbe, lw=3 if k == n_sel else 1.5)
        ax.text(1.03, E_k, f"n={k}", va="center", fontsize=9, color=farbe)
    ax.set_ylabel("Energie [eV]")
    ax.set_xticks([])
    ax.set_xlim(0, 1.3)
    ax.spines[["top", "right", "bottom"]].set_visible(False)
    st.pyplot(fig)

    st.markdown(
        """
**Das ist eine Treppe ohne Rampe.** Zwischen den Stufen ist nichts. Nicht
„schwer erreichbar", sondern schlicht nicht existent. Das Teilchen kann diese
Energien haben und keine anderen.

Genau dieses Wort, dass Energie in Stufen kommt statt fließend, ist der Grund
für den Namen: **Quant** heißt auf Latein „wie viel", ein Quantum ist eine
abgezählte Portion.
"""
    )

    st.divider()
    st.subheader("Zwei Dinge, die man sofort sieht")

    st.markdown(
        """
**Erstens: Die Stufen werden nach oben immer weiter auseinander.**
Die Energie wächst mit n². Von n=1 auf n=2 vervierfacht sie sich, von n=1 auf
n=3 verneunfacht sie sich.

**Zweitens: Enger Kasten = größere Abstände.** Zieh den Breiten-Regler klein und
schau, wie die Leiter auseinanderzieht.
"""
    )
    st.warning(
        "**Das ist der Satz, auf den es ankommt:**\n\n"
        "Je enger du ein Teilchen einsperrst, desto heftiger wehrt es sich. "
        "Ein eingesperrtes Teilchen kann nicht einfach stillstehen, es hat "
        "immer eine Mindestenergie, die sogenannte Nullpunktsenergie.\n\n"
        "Und genau deshalb stürzen Elektronen nicht in den Atomkern. "
        "Ein Elektron im Kern wäre extrem eng eingesperrt und bräuchte dafür "
        "absurd viel Energie. Die hat es nicht. Also bleibt es draußen."
    )


# ==================================================================
# Kapitel 4
# ==================================================================
elif kapitel == KAPITEL[4]:
    st.title("Farbe aus dem Sprung")
    st.markdown(
        """
Jetzt wird es sichtbar. Wenn ein Teilchen von einer Stufe auf eine tiefere
fällt, muss die Energiedifferenz irgendwo hin. Sie wird als **Lichtteilchen**
(Photon) abgegeben.

Und die Energie des Photons bestimmt seine **Farbe**. Große Differenz = blau,
kleine Differenz = rot.

Das heißt: Aus der Treppe wird Licht, das man messen kann. **So verifiziert man
diese ganze Theorie überhaupt.**
"""
    )

    c1, c2, c3 = st.columns(3)
    with c1:
        n_hoch = st.selectbox("Von n =", [2, 3, 4, 5, 6], index=0)
    with c2:
        n_tief = st.selectbox("Nach n =", [1, 2, 3, 4, 5], index=0)
    with c3:
        L_nm = st.slider("L [nm]", 0.4, 2.0, 0.7, 0.05, key="l4")

    if n_hoch <= n_tief:
        st.error("Der Startzustand muss höher liegen als der Zielzustand.")
    else:
        L = L_nm * 1e-9
        dE = energie(n_hoch, L) - energie(n_tief, L)
        wl_nm = (H * C / dE) * 1e9
        farbe = wellenlaenge_zu_farbe(wl_nm)

        m1, m2 = st.columns(2)
        m1.metric("Freigesetzte Energie", f"{dE / EV:.2f} eV")
        m2.metric("Wellenlänge", f"{wl_nm:.0f} nm")

        fig, ax = plt.subplots(figsize=(7, 1.2))
        ax.add_patch(plt.Rectangle((0, 0), 1, 1, color=farbe))
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
        st.pyplot(fig)

        if 380 <= wl_nm <= 750:
            st.success(f"Sichtbares Licht bei {wl_nm:.0f} nm.")
        elif wl_nm < 380:
            st.info(f"{wl_nm:.0f} nm – Ultraviolett, für dein Auge unsichtbar.")
        else:
            st.info(f"{wl_nm:.0f} nm – Infrarot, spürbar als Wärme.")

    st.divider()
    st.subheader("Warum das kein Spielzeug ist")
    st.markdown(
        """
Spiel mit dem Breiten-Regler und beobachte, wie sich die Farbe ändert.

Winzige Kristalle von wenigen Nanometern Größe heißen Quantenpunkte. Ihre Farbe
hängt nicht von ihrem Material ab, sondern von ihrer **Größe**, genau nach dem
Prinzip, das du hier gerade verstellst.

Kleiner Punkt → enger Kasten → größere Sprünge → blaueres Licht.

In QLED-Fernsehern steckt genau das. Die Kastenbreite, an der du gerade drehst,
ist dort eine Fertigungstoleranz in der Produktion.
"""
    )
    st.info(
        "**Für den Diskursraum:** Zwischen Schrödingers Gleichung (1926) und dem "
        "Fernseher im Wohnzimmer liegen rund 90 Jahre. Bei Haber-Bosch lagen "
        "zwischen Laborreaktion und Weltproduktion keine zehn. "
        "Warum eigentlich so unterschiedlich?"
    )


# ==================================================================
# Kapitel 5
# ==================================================================
elif kapitel == KAPITEL[5]:
    st.title("Wo ist das Teilchen?")
    st.markdown(
        """
Hier kommt der Bruch mit dem alten Weltbild. Bis jetzt klang alles noch
harmlos – Wellen, Stufen, Licht. Jetzt wird es unangenehm.

**Frage: Wo genau ist das Elektron im Kasten?**

Antwort der klassischen Physik: An einem bestimmten Ort. Wir wissen ihn
vielleicht nicht, aber es gibt ihn.

Antwort der Quantenmechanik: **Die Frage hat keine Antwort.** Es gibt keinen
Ort, an dem es ist. Es gibt nur die Wahrscheinlichkeit, es irgendwo zu finden,
wenn man nachschaut.
"""
    )

    c1, c2 = st.columns(2)
    with c1:
        n = st.slider("Quantenzahl n", 1, 6, 2, key="n5")
    with c2:
        modus = st.radio(
            "Ansicht",
            ["Welle (ψ)", "Wahrscheinlichkeit (ψ²)"],
            horizontal=True,
        )

    st.pyplot(
        plot_zustand(
            n, 1.0, "wahrscheinlichkeit" if "Wahrsch" in modus else "welle"
        )
    )

    if n >= 2 and "Wahrsch" in modus:
        st.warning(
            f"**Schau dir die Nullstellen an.** Bei n={n} gibt es {n-1} "
            f"{'Stelle' if n == 2 else 'Stellen'} im Kasten, an denen das "
            "Teilchen **nie** gefunden wird.\n\n"
            "Wie kommt das Teilchen nun von links nach rechts, wenn es in der Mitte niemals sein "
            "kann?\n\n"
            "Die Antwort ist nicht „es springt darüber\". Die Antwort ist, dass "
            "die Frage eine falsche Vorstellung enthält, nämlich die, dass das "
            "Teilchen eine Bahn hat, auf der es entlangläuft. Hat es aber nicht."
        )

    st.divider()
    st.markdown(
        """
### Was hier wirklich passiert ist

Die Quantenmechanik nimmt einen Begriff aus dem Weltbild heraus, den vorher
niemand für verhandelbar gehalten hätte: **dass Dinge einen Ort haben.**

Das ist keine Messungenauigkeit. Es ist keine Frage besserer Geräte. Nach
allem, was wir wissen, hat das Elektron zwischen zwei Messungen schlicht keinen
definierten Ort.

Einstein hat das bis zu seinem Tod nicht akzeptiert. Sein berühmter Einwand,
Gott würfle nicht, war kein Scherz, sondern ein ernst gemeinter Protest gegen
genau diesen Punkt. Die Experimente haben ihm später unrecht gegeben.
"""
    )
    st.info(
        "**Die Transformation:** Hier verändert sich nicht das "
        "Wissen über die Welt, sondern die Vorstellung davon, welche Fragen "
        "überhaupt sinnvoll sind. Das ist eine andere Art von Umbruch als bei "
        "Haber-Bosch."
    )


# ==================================================================
# Kapitel 6
# ==================================================================
elif kapitel == KAPITEL[6]:
    st.title("Warum merkst du davon nichts?")
    st.markdown(
        """
Berechtigter Einwand: warum ist die Welt um dich herum dann so normal?

Machen wir die Probe. Wir sperren statt eines Elektrons **Dich** in einen Raum.
Dieselbe Formel, andere Zahlen.
"""
    )

    c1, c2 = st.columns(2)
    with c1:
        masse = st.number_input("Dein Gewicht [kg]", 30.0, 200.0, 70.0, 5.0)
    with c2:
        raum = st.number_input("Raumbreite [m]", 1.0, 20.0, 5.0, 0.5)

    E1 = energie(1, raum, masse)
    E2 = energie(2, raum, masse)
    sprung = E2 - E1

    st.metric("Dein Sprung von n=1 auf n=2", f"{sprung:.2e} Joule")

    st.markdown(
        f"""
Zum Vergleich: Ein Sandkorn einen Millimeter anzuheben kostet ungefähr
0,00000001 Joule. Das ist etwa **{1e-8 / sprung:.0e} mal mehr** als dein
Quantensprung.

Deine Energiestufen sind also da. Sie sind nur so absurd eng beieinander, dass
keine Messung der Welt sie unterscheiden könnte. Die Treppe ist noch da,
aber die Stufen sind niedriger als die Rauheit des Materials.
"""
    )
    st.success(
        "**Das ist ein wichtiger Punkt gegen ein verbreitetes Missverständnis:**\n\n"
        "Die Quantenmechanik hat die klassische Physik nicht widerlegt. "
        "Sie hat sie **eingeordnet**. Newton ist nicht falsch, sondern ein "
        "Grenzfall für große, schwere Dinge. Und der stimmt dort so gut, dass "
        "wir bis heute Raumsonden damit steuern.\n\n"
        "Wissenschaftlicher Fortschritt heißt hier nicht: das Alte war Unsinn. "
        "Sondern: wir wissen jetzt, wo seine Grenzen liegen."
    )
    st.markdown(
        """
Eine Transformation des Weltbilds bedeutet selten, dass alles Vorherige weg ist.
Meistens wird das Alte zu einem **Spezialfall** des Neuen.
"""
    )


# ==================================================================
# Kapitel 7
# ==================================================================
else:
    st.title("Diskussion")
    st.markdown(
        "Vier Fragen für die interdisziplinäre Runde. "
    )

    fragen = [
        ("Was passiert mit einem Weltbild, wenn es kippt?",
         "1898 hätte jeder Physiker gesagt: Energie ist eine kontinuierliche "
         "Größe, und Dinge haben einen Ort. Dreißig Jahre später war beides "
         "nicht mehr haltbar.\n\n"
         "Woher wollen wir wissen, dass unsere heutigen "
         "Selbstverständlichkeiten stabiler sind? Und was folgt daraus für "
         "Fächer, die keine Experimente machen können?"),
        ("Erkenntnis oder Anwendung, was transformiert eigentlich?",
         "Haber-Bosch veränderte binnen weniger Jahre, wie viele Menschen die "
         "Erde ernähren kann und zugleich, wie Kriege geführt werden. Die "
         "Quantenmechanik veränderte zunächst nur, wie wir die Welt denken.\n\n"
         "Ist eine Erkenntnis, die nichts verändert, weniger transformativ? "
         "Oder mehr, weil sie tiefer sitzt?"),
        ("Wer trägt Verantwortung für die Folgen?",
         "Schrödinger und Bohr haben keine Fernseher, keine Solarzellen und "
         "keine Atombomben gebaut. Ihre Gleichungen stehen aber hinter allen "
         "dreien.\n\n"
         "Endet Verantwortung an der Labortür? Und falls nicht, wie weit "
         "reicht sie, wenn zwischen Erkenntnis und Anwendung Jahrzehnte "
         "liegen?"),
        ("Was heißt Verstehen, wenn Anschauung versagt?",
         "Wir haben heute Bilder benutzt: Saiten, Kästen, Treppen. Alle sind "
         "streng genommen falsch. Richtig ist nur die Mathematik.\n\n"
         "Ist das ein Problem? Kann man etwas verstehen, das man sich nicht "
         "vorstellen kann? Und gibt es das in euren Fächern auch?"),
    ]

    for titel, text in fragen:
        with st.expander(titel):
            st.write(text)

    st.divider()
    st.caption(
        "Arbeitskreis „Vom Dünger zur Quantenwelt – wie Wissenschaft unser "
        "Weltbild transformiert“ · Constantin Richard Feitl & Dato Tsomaia"
    )
