# -*- coding: utf-8 -*-
"""
Vom Dünger zur Quantenwelt – Arbeitskreis
Block: Haber-Bosch-Verfahren (Quantenchemie-Fokus)
 
Start:   python -m streamlit run haber_bosch.py
Braucht: streamlit, numpy, matplotlib
         sowie hf_pure.py im selben Ordner (reines numpy, kein Compiler
         noetig - laeuft nativ unter Windows, macOS, Linux)
 
    pip install streamlit numpy matplotlib
 
"""
 
import io
import contextlib
import traceback
 
import numpy as np
import matplotlib.pyplot as plt
import streamlit as st
 
try:
    import hf_pure
    HF_DA = True
except ImportError:
    HF_DA = False
 
# ==================================================================
# Konstanten
# ==================================================================
BLAU = "#1b6ca8"
ORANGE = "#d95f02"
GRUEN = "#2a9d5c"
GRAU = "#c9c9c9"
 
HARTREE_KJ = 2625.499639  # 1 Hartree in kJ/mol
 
st.set_page_config(page_title="Haber-Bosch & Quantenchemie", layout="centered")
 
 
# ==================================================================
# Quantenchemie-Hilfsfunktionen (eigene Hartree-Fock-Implementierung,
# siehe hf_pure.py - reines numpy, STO-3G, Elemente H, C, N, O)
# ==================================================================
@st.cache_data(show_spinner=False)
def qm_energie(atome, spin=0):
    """Gesamtenergie eines Moleküls per Hartree-Fock (STO-3G).
 
    atome: Tupel aus (Elementsymbol, (x,y,z))-Paaren, Angstrom, z.B.
           (("H", (0,0,0)), ("H", (0,0,0.74)))
    spin:  Anzahl ungepaarter Elektronen (0 für die meisten Moleküle
           hier, 1 für ein einzelnes H-Atom)
    Gibt die Energie in Hartree zurück.
    """
    return hf_pure.energie(list(atome), spin=spin)
 
 
# Pyramidale NH3-Geometrie (Angstrom), N-H = 1.012 A, Winkel = 106.7°
NH3_GEOMETRIE = (
    ("N", (0.000000, 0.000000, 0.116489)),
    ("H", (0.000000, 0.939731, -0.271808)),
    ("H", (0.813831, -0.469865, -0.271808)),
    ("H", (-0.813831, -0.469865, -0.271808)),
)
 
# Vorberechnete Bindungskurven (Hartree-Fock/6-31G) fuer den Teaser in
# Kapitel 1 - damit die App dort ohne Wartezeit startet.
SCAN_N2_R = [0.90, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60,
             1.70, 1.80, 1.90, 2.00, 2.10, 2.20]
SCAN_N2_E = [494.3, 85.0, 0.0, 83.6, 247.2, 441.1, 639.3, 829.8,
             1008.0, 1173.1, 1325.3, 1465.1, 1593.0, 1709.6]
SCAN_H2_R = [0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30,
             1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.10, 2.20]
SCAN_H2_E = [178.8, 42.3, 0.0, 6.3, 37.9, 82.2, 132.4, 184.7, 237.1,
             288.3, 337.7, 385.0, 430.1, 472.8, 513.1, 551.0, 586.5, 619.8]
 
 
# ==================================================================
# Der Code-Runner: hier wird selbst gerechnet
# ==================================================================
def rechne_selbst(key, code, hinweis=None, hoehe=260):
    """Editierbares Codefeld mit Ausführen-Knopf und Ausgabe."""
    st.markdown("##### Rechnung selbst ausführen")
    if hinweis:
        st.caption(hinweis)
 
    if f"code_{key}" not in st.session_state:
        st.session_state[f"code_{key}"] = code
 
    eingabe = st.text_area(
        "Python-Code – du darfst alles ändern",
        value=st.session_state[f"code_{key}"],
        height=hoehe,
        key=f"area_{key}",
        label_visibility="collapsed",
    )
 
    c1, c2 = st.columns([1, 3])
    with c1:
        los = st.button("▶ Ausführen", key=f"run_{key}", type="primary")
    with c2:
        if st.button("Zurücksetzen", key=f"reset_{key}"):
            st.session_state[f"code_{key}"] = code
            st.rerun()
 
    if los:
        puffer = io.StringIO()
        umgebung = {"np": np, "numpy": np, "plt": plt}
        if HF_DA:
            umgebung.update({
                "qm_energie": lambda atome, spin=0: hf_pure.energie(list(atome), spin=spin),
                "HARTREE_KJ": HARTREE_KJ,
                "NH3_GEOMETRIE": NH3_GEOMETRIE,
            })
        try:
            with contextlib.redirect_stdout(puffer):
                exec(eingabe, umgebung)
            ausgabe = puffer.getvalue()
            if ausgabe.strip():
                st.code(ausgabe, language=None)
            else:
                st.info("Kein print() im Code – deshalb keine Ausgabe.")
            fig = plt.gcf()
            if fig.get_axes():
                st.pyplot(fig)
            plt.close("all")
        except Exception:
            st.error("Da ist etwas schiefgegangen:")
            st.code(traceback.format_exc().splitlines()[-1], language=None)
 
 
# ==================================================================
# Navigation
# ==================================================================
KAPITEL = [
    "0 · Wir stehen in einem Meer aus Dünger",
    "1 · Warum Luft nicht düngt",
    "2 · Moleküle aus dem Nichts berechnen",
    "3 · Brot und Sprengstoff",
    "4 · Diskussion",
]
 
st.sidebar.title("Vom Dünger zur Quantenwelt")
st.sidebar.caption("Block Haber-Bosch")
kapitel = st.sidebar.radio("Kapitel", KAPITEL, label_visibility="collapsed")
st.sidebar.divider()
st.sidebar.markdown(
    "**Zum Mitrechnen:** In manchen Kapiteln steht ein Codefeld. "
    "Du kannst Zahlen ändern und auf *Ausführen* drücken. "
    "Kaputtmachen kannst du nichts, *Zurücksetzen* holt das Original zurück."
)
if not HF_DA:
    st.sidebar.error(
        "`hf_pure.py` wurde nicht gefunden. Kapitel 2 braucht diese Datei "
        "im selben Ordner wie `haber_bosch.py`."
    )
 
 
# ==================================================================
# Kapitel 0
# ==================================================================
if kapitel == KAPITEL[0]:
    st.title("Wir stehen in einem Meer aus Dünger")
    st.markdown(
        """
Die Luft in diesem Raum besteht zu **78 Prozent aus Stickstoff**. Mit jedem
Atemzug ziehst du Milliarden Stickstoffmoleküle ein und atmest sie unverändert
wieder aus.
 
Gleichzeitig ist Stickstoff der Nährstoff, an dem Pflanzenwachstum als Erstes
scheitert. Er steckt in jedem Eiweiß, in jedem Stück DNA, in jedem Muskel.
 
Das ist die Ausgangslage, und sie ist absurd: **Der Rohstoff ist überall, und
niemand kommt an ihn heran.**
"""
    )
    st.info(
        "Über Jahrtausende war Stickstoff deshalb ein Engpass. Man sammelte "
        "Mist, holte Guano von südamerikanischen Inseln, baute Salpeter in der "
        "Atacama-Wüste ab. Um 1900 war absehbar, dass diese Quellen für eine "
        "wachsende Weltbevölkerung nicht reichen würden."
    )
 
    m1, m3 = st.columns(2)
    m1.metric("Stickstoff in der Luft", "78 %")
    m3.metric("Menschen heute davon ernährt", "≈ 50 %")
 
    st.divider()
    st.markdown(
        """
1909 löste Fritz Haber das Problem im Labor, 1913 baute Carl Bosch daraus eine
Fabrik. Die Reaktion, um die sich alles dreht, ist auf den ersten Blick simpel:
"""
    )
    st.latex(r"\mathrm{N_2} + 3\,\mathrm{H_2} \;\longrightarrow\; 2\,\mathrm{NH_3}")
    st.markdown(
        "Warum das trotzdem so schwer war, liegt an etwas, das man nicht "
        "sieht: wie fest die beiden Stickstoffatome in N₂ aneinanderhängen. "
        "Darum geht es im nächsten Kapitel – und danach rechnen wir es selbst."
    )
 
    st.divider()
    rechne_selbst(
        "luft",
        '''# Wie viel Stickstoff ist in diesem Raum?
# Aendere die Raumgroesse und schau, was rauskommt.
 
laenge = 10     # Meter
breite = 8
hoehe  = 3
 
volumen = laenge * breite * hoehe        # Kubikmeter
luft_kg = volumen * 1.2                  # Luft wiegt 1,2 kg pro Kubikmeter
stickstoff_kg = luft_kg * 0.755          # 75,5 Massenprozent
 
print(f"Raumvolumen:        {volumen} m3")
print(f"Luft im Raum:       {luft_kg:.0f} kg")
print(f"Davon Stickstoff:   {stickstoff_kg:.0f} kg")
print()
 
# Ein Hektar Weizen braucht ungefaehr 150 kg Stickstoff pro Jahr.
print(f"Das reicht rechnerisch fuer {stickstoff_kg/150:.1f} Hektar Weizen.")
print("Wenn man drankaeme. Genau das ist das Problem.")
''',
        hinweis="Trag die ungefähre Größe eures Seminarraums ein.",
        hoehe=300,
    )
 
 
# ==================================================================
# Kapitel 1
# ==================================================================
elif kapitel == KAPITEL[1]:
    st.title("Warum Luft nicht düngt")
    st.markdown(
        """
Stickstoff kommt in der Luft nie einzeln vor, sondern immer paarweise: N₂.
Die beiden Atome halten sich mit einer **Dreifachbindung** fest, die sehr stabil ist.
 
Die Bindungsenergie von N≡N ist mehr als doppelt so groß wie die Bindungsenergie von H-H:
"""
    )
 
    v1, v2 = st.columns(2)
    v1.metric("Bindungsenergie N≡N", "945 kJ/mol")
    v2.metric("Bindungsenergie H–H", "436 kJ/mol")
 
    st.warning(
        "**Das ist der ganze Grund für die Schwierigkeit.** N₂ ist chemisch "
        "fast so träge wie ein Edelgas. Es reagiert bei Raumtemperatur "
        "praktisch mit nichts."
    )
 
 
# ==================================================================
# Kapitel 2 – Quantenchemie live
# ==================================================================
elif kapitel == KAPITEL[2]:
    st.title("Moleküle aus dem Nichts berechnen")
    st.markdown(
        """
 
Man kann Bindungsenergien auch quantenmechanisch ausrechnen. Ganz ohne Experiment. Man
gibt nur an: Wo stehen die Atomkerne, wie viele Elektronen gibt es. Alles
Weitere folgt aus der Schrödinger-Gleichung.
 
Das ist ein bemerkenswerter Anspruch: **Chemie als Rechenaufgabe.**
"""
    )
 
    if not HF_DA:
        st.error(
            "Die Datei `hf_pure.py` liegt nicht neben dieser App. Leg sie "
            "in denselben Ordner, dann funktioniert dieses Kapitel."
        )
        st.stop()
 
    with st.expander("**Wie rechnet der Computer das aus?** (eine Analogie)"):
        st.markdown(
            """
Ein Molekül ist ein Haufen Elektronen, die sich alle gegenseitig abstoßen.
Exakt lösen lässt sich das nicht mehr, sobald es mehr als zwei sind.
 
Der Trick, den wir hier benutzen (**Hartree-Fock**): Stell dir eine volle
Party vor. Statt zu berechnen, wie jeder Gast jedem einzelnen anderen
ausweicht, tut man so, als bewege sich jeder nur durch eine **gleichmäßige
Menschenmenge**. Das ist eine Näherung: grob, aber schnell. Aus Demonstrationszwecken genügt hier Hartree-Fock.
"""
        )
 
    st.divider()
    st.subheader("Schritt 1 · Ein einzelner Punkt")
    st.markdown(
        "Wir nehmen das einfachste Molekül überhaupt: **H₂**, zwei "
        "Wasserstoffatome. Wir geben nur vor, wie weit die Kerne "
        "auseinanderstehen, der Computer liefert die Energie."
    )
 
    abstand = st.slider("Abstand der beiden H-Kerne [Ångström]",
                        0.40, 2.20, 0.74, 0.02)
 
    if st.button("Energie berechnen", type="primary"):
        with st.spinner("Schrödinger-Gleichung wird gelöst …"):
            E = qm_energie((("H", (0, 0, 0)), ("H", (0, 0, abstand))))
        st.session_state["h2_punkt"] = E
 
    if "h2_punkt" in st.session_state:
        E = st.session_state["h2_punkt"]
        c1, c2 = st.columns(2)
        c1.metric("Gesamtenergie", f"{E:.4f} Hartree")
        c2.metric("In vertrauteren Einheiten", f"{E * HARTREE_KJ:,.0f} kJ/mol")
        st.caption(
            "„Hartree“ ist die natürliche Energieeinheit der Quantenchemie."
        )
 
    st.divider()
    st.subheader("Schritt 2 · Die ganze Bindungskurve")
    st.markdown(
        "Jetzt wiederholen wir das für 20 verschiedene Abstände – dieselbe "
        "Kurve wie in Kapitel 1, nur diesmal **live selbst gerechnet** statt "
        "vorgegeben."
    )
 
    if st.button("Kurve berechnen (ein paar Sekunden)"):
        with st.spinner("Die Schrödinger-Gleichung wird 20 Mal gelöst …"):
            rs = np.linspace(0.4, 2.2, 20)
            es = [qm_energie((("H", (0, 0, 0)), ("H", (0, 0, float(r)))))
                  for r in rs]
            E_atom = qm_energie((("H", (0, 0, 0)),), spin=1)
        st.session_state["h2_kurve"] = (rs, es, E_atom)
 
    if "h2_kurve" in st.session_state:
        rs, es, E_atom = st.session_state["h2_kurve"]
        i = int(np.argmin(es))
        fig, ax = plt.subplots(figsize=(7, 3.4))
        ax.plot(rs, (np.array(es) - min(es)) * HARTREE_KJ, "o-",
                color=BLAU, lw=2, ms=4)
        ax.scatter([rs[i]], [0], color=ORANGE, s=90, zorder=5)
        ax.axvline(0.741, color=GRUEN, ls="--", lw=1.5)
        ax.text(0.76, 20, "gemessen: 0.741 Å", color=GRUEN, fontsize=9)
        ax.set_xlabel("Kernabstand [Ångström]")
        ax.set_ylabel("Energie relativ zum Minimum [kJ/mol]")
        ax.spines[["top", "right"]].set_visible(False)
        st.pyplot(fig)
 
        bindung = (2 * E_atom - es[i]) * HARTREE_KJ
        d1, d2 = st.columns(2)
        d1.metric("Berechnete Bindungslänge", f"{rs[i]:.2f} Å",
                  delta=f"{rs[i] - 0.741:+.2f} Å gegen Messung")
        d2.metric("Berechnete Bindungsstärke", f"{bindung:.0f} kJ/mol",
                  delta=f"{bindung - 436:+.0f} gegen Messung")
 
        st.success(
            "**Kurz innehalten.** Wir haben nichts gemessen. Wir haben nur "
            "gesagt: hier sind zwei Protonen und zwei Elektronen. Heraus kam, "
            "wie lang die Bindung ist und wie stark sie hält, beides nah an "
            "den echten Werten aus dem Labor."
        )
 
    st.divider()
    st.subheader("Schritt 3 · Jetzt die eigentliche Haber-Bosch-Reaktion")
    st.markdown(
        "Derselbe Trick, dreimal angewandt – auf N₂, H₂ und NH₃ und wir "
        "können ausrechnen, ob die Reaktion Energie freisetzt. Genau die "
        "Frage, die am Anfang stand."
    )
    st.latex(r"\mathrm{N_2} + 3\,\mathrm{H_2} \;\longrightarrow\; 2\,\mathrm{NH_3}")
 
    st.caption(
        "Rechnung mit STO-3G, dem kleinsten gebräuchlichen Basissatz. "
        "Größere Basissätze wären genauer, aber auch deutlich rechenaufwendiger."
    )
 
    if st.button("Reaktionsenergie berechnen", type="primary"):
        with st.spinner("Drei Moleküle werden durchgerechnet (ein paar Sekunden) …"):
            E_N2 = qm_energie((("N", (0, 0, 0)), ("N", (0, 0, 1.10))))
            E_H2 = qm_energie((("H", (0, 0, 0)), ("H", (0, 0, 0.74))))
            E_NH3 = qm_energie(NH3_GEOMETRIE)
            dE = (2 * E_NH3 - (E_N2 + 3 * E_H2)) * HARTREE_KJ
        st.session_state["hb_ergebnis"] = dE
 
    if "hb_ergebnis" in st.session_state:
        dE = st.session_state["hb_ergebnis"]
        st.metric(
            "Berechnete Reaktionsenergie",
            f"{dE:+.0f} kJ/mol",
            delta=f"{dE - (-92.4):+.0f} kJ/mol gegen Messung (−92,4)",
            delta_color="inverse",
        )
 
        if dE < 0:
            st.success(
                "**Das Vorzeichen stimmt: Energie wird frei.** Die Reaktion "
                "läuft also grundsätzlich in die richtige Richtung, genau "
                "wie Haber es 1909 im Labor gefunden hat, nur dass wir es "
                "gerade eben am Computer nachvollzogen haben."
            )
        st.info(
            "**Der genaue Zahlenwert trifft die Messung "
            "nicht exakt. STO-3G ist der kleinste gebräuchliche Basissatz."
            "eine sehr grobe Beschreibung der Elektronenwolken. Größere "
            "Basissätze und bessere Methoden kämen näher an die −92,4 "
            "kJ/mol heran, kosten aber massiv mehr Rechenzeit. **Das "
            "Vorzeichen und die Größenordnung stimmen trotzdem** – und "
            "genau das haben wir gerade selbst nachgerechnet, ohne Labor."
        )
 
    st.divider()
    rechne_selbst(
        "qm",
        '''# Spiel selbst mit der Quantenchemie.
# qm_energie(atome, spin=0) gibt die Energie in Hartree zurueck.
# atome ist eine Liste aus (Elementsymbol, (x,y,z)) in Angstrom.
# Verfuegbare Elemente: H, C, N, O.
# 1 Hartree = 2625.5 kJ/mol (steht schon in HARTREE_KJ).
 
# Beispiel: Bindungslaenge von CO (Kohlenmonoxid) grob suchen
for abstand in [0.9, 1.0, 1.1, 1.128, 1.2, 1.3]:
    geometrie = [("C", (0, 0, 0)), ("O", (0, 0, abstand))]
    E = qm_energie(geometrie)
    print(f"R = {abstand:.3f} A   ->   E = {E:+.4f} Hartree")
 
print()
print("Das Minimum zeigt die berechnete Bindungslaenge.")
print("Gemessener Wert fuer C-O: 1.128 Angstrom.")
''',
        hinweis="Ändere Elemente (H, C, N, O), Positionen oder Abstände.",
        hoehe=340,
    )
 
 
# ==================================================================
# Kapitel 3
# ==================================================================
elif kapitel == KAPITEL[3]:
    st.title("Brot und Sprengstoff")
    st.markdown(
        """
Bis hierher war es eine Erfolgsgeschichte: ein Engpass, eine Reaktion, die
sich am Computer nachrechnen lässt, eine Fabrik. Jetzt der andere Teil.
 
**Ammoniak ist der Ausgangsstoff für beides.** Aus NH₃ wird über Salpetersäure
sowohl Kunstdünger als auch Sprengstoff. Derselbe Prozess, dieselbe Anlage,
dieselben Menschen.
"""
    )
 
    st.divider()
    rechne_selbst(
        "menschen",
        '''# Wie viele Menschen haengen an diesem Verfahren?
 
weltbevoelkerung = 8.1        # Milliarden
anteil_kunstduenger = 0.50    # Schaetzung nach Vaclav Smil
 
menschen = weltbevoelkerung * anteil_kunstduenger
print(f"Rund {menschen:.1f} Milliarden Menschen essen Nahrung,")
print(f"deren Stickstoff aus diesem Verfahren stammt.")
print()
 
# Und die andere Seite:
ammoniak_welt = 180           # Millionen Tonnen pro Jahr
anteil_duenger = 0.80
 
print(f"Weltproduktion Ammoniak:  {ammoniak_welt} Mio. Tonnen/Jahr")
print(f"davon fuer Duenger:       {ammoniak_welt*anteil_duenger:.0f} Mio. t")
print(f"fuer alles andere:        {ammoniak_welt*(1-anteil_duenger):.0f} Mio. t")
print()
print("'Alles andere' heisst: Kunststoffe, Reinigungsmittel - und Sprengstoff.")
''',
        hinweis="Verschieb den Anteil und überleg: Ab welchem Wert wird die "
                "Abschaffung dieses Verfahrens undenkbar?",
        hoehe=400,
    )
 
    st.divider()
    st.markdown(
        """
### Die andere Seite
 
Vor 1914 importierte das Deutsche Reich seinen Salpeter aus Chile. Mit der
britischen Seeblockade war dieser Weg zu. Nach gängiger Einschätzung hätte der
Munitionsnachschub binnen etwa eines Jahres geendet.
 
Die Ammoniakanlagen der BASF ersetzten diesen Import. **Der Erste Weltkrieg
konnte auch deshalb jahrelang weitergeführt werden** (vgl. Szöllösi-Janze 1998,
S. 270f.).
 
Fritz Haber ging darüber hinaus. Er organisierte den Einsatz chemischer
Kampfstoffe an der Front und war bei Ypern 1915 persönlich anwesend
(vgl. Szöllösi-Janze 1998, S. 320f.). Seine Frau Clara Immerwahr, selbst
promovierte Chemikerin, nahm sich wenige Tage danach das Leben.
 
1918 erhielt Haber den Nobelpreis für Chemie – für die Ammoniaksynthese.
1933 musste er als Jude aus Deutschland emigrieren. Er starb 1934 im Exil.
"""
    )
    st.warning(
        "**Die Versuchung ist, sich für eine Seite zu entscheiden.** Held oder "
        "Kriegsverbrecher, Ernährer oder Giftgasorganisator. Beides ist "
        "belegbar, und beides ist zu einfach.\n\n"
        "Interessanter ist die Frage, was der Fall über die Struktur "
        "wissenschaftlicher Arbeit sagt – und ob eine andere Person an Habers "
        "Stelle anders gehandelt hätte."
    )
 
 
# ==================================================================
# Kapitel 4
# ==================================================================
else:
    st.title("Diskussion")
    st.markdown("Vier Fragen.")
 
    fragen = [
        ("Ist die Erkenntnis schuldig oder erst die Anwendung?",
         "Die Reaktionsgleichung N₂ + 3 H₂ → 2 NH₃ ist wertfrei. Aus dem "
         "Produkt wird Brot oder Munition, je nachdem, wer die Anlage "
         "besitzt.\n\n"
         "Lässt sich diese Trennung durchhalten? Oder ist sie eine bequeme "
         "Ausrede, zumal Haber die militärische Verwendung nicht nur hinnahm, "
         "sondern selbst betrieb?"),
        ("Hätte man es lassen können?",
         "Angenommen, Haber hätte 1909 abgebrochen. Das Problem wurde weltweit "
         "bearbeitet, die Vorarbeiten lagen vor, der Bedarf war enorm.\n\n"
         "Wenn eine Entdeckung ohnehin fällig ist, ändert individueller "
         "Verzicht dann überhaupt etwas? Und falls nicht: bleibt trotzdem eine "
         "persönliche Verantwortung?"),
        ("Was ist mit den Folgen, die niemand wollte?",
         "Der Stickstoff aus diesen Anlagen landet heute in Grundwasser, in "
         "überdüngten Flüssen und als Lachgas in der Atmosphäre, ein "
         "Treibhausgas. Diese Folgen hat weder Haber noch Bosch beabsichtigt "
         "oder auch nur gekannt.\n\n"
         "Kann man für etwas verantwortlich sein, das zum Zeitpunkt der "
         "Handlung nicht absehbar war? Hans Jonas hat genau daraus eine neue "
         "Ethik abgeleitet."),
        ("Warum so schnell – und die Quantenmechanik so langsam?",
         "Zwischen Habers Laborreaktion (1909) und der industriellen Anlage "
         "(1913) lagen vier Jahre. Zwischen Schrödingers Gleichung (1926) und "
         "den Rechnungen, die wir gerade eben in Sekunden gemacht haben, "
         "lagen Jahrzehnte Entwicklung an Theorie und Rechenleistung.\n\n"
         "Woran liegt dieser Unterschied? An der Wissenschaft, an der "
         "wirtschaftlichen Nachfrage, am Krieg? Und was heißt das für die "
         "Frage, wann Forschung überhaupt steuerbar ist?"),
    ]
 
    for titel, text in fragen:
        with st.expander(titel):
            st.write(text)
 
    st.divider()
    st.caption(
        "Arbeitskreis „Vom Dünger zur Quantenwelt – wie Wissenschaft unser "
        "Weltbild transformiert“ · Constantin Richard Feitl & Dato Tsomaia"
    )
