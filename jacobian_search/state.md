# Jacobian Search State - Részletes Napló

## Cél
Keresni F1 = x + A, F2 = y + B, A,B ∈ (x,y)^2, det J ≡1, A(1,0)=-1, B(1,0)=0.
J feltétel: A_x + B_y + A_x B_y - A_y B_x =0

## Jelölt családok és motiváció

### 1. Teljes együttható-illesztés d=3,4,5
Motiváció: Brute force minden együtthatóra, Gröbner-bázis, elimináció.
- d=3: 14 ismeretlen, 14+2 egyenlet. Homogén dekompozícióval redukálva 9 ismeretlen, 11 egyenlet.
- d=4: 24 ismeretlen, 27+2 egyenlet, redukálva 15 ismeretlen, 20 egyenlet.
- d=5: 36 ismeretlen, 44+2 egyenlet, redukálva 22 ismeretlen, 32 egyenlet.

### 2. Ritka támogatások szisztematikusan változtatott monom-mintákkal
Motiváció: Kevés monom -> kevesebb ismeretlen, de sok egyenlet, overdetermined, de kis Gröbner könnyebb.
- d=3: monomok 7/db, support méret 2: 21 db, párok 441, 200 kipróbálva, 0 megoldás.
- d=3 support méret 3: 35 db, random 100, 0 megoldás.
- d=4 support méret 4: 495 db, random 100, 0 megoldás.
- Következtetés: ritka támogatások nem adnak megoldást d=3,4-re.

### 3. Súlyozott fokú és vegyes fokú konstrukciók
Motiváció: Súlyozott homogén Keller térképek osztályozása.
- Lemma: Ha A,B súlyozott homogének azonos súlyozott fokkal és w_x != w_y, akkor A_x és B_y különböző súlyozott fokúak, így külön-külön nullának kell lenniük -> A=A(y), B=B(x) -> automorfizmus, nincs ütközés.
- Következtetés: Súlyozott homogén család nem ad ellenpéldát.

### 4. Strukturált polinomkompozíciók és koordinátatranszformációk
Motiváció: Drużkowski térképek: F_i = x_i + (lineáris forma)^3, JH nilpotens.
- 2D Drużkowski: F1 = x + (a x + b y)^3, F2 = y + (c x + d y)^3
- Feltételek: a^3 + d c^2=0, a^2 b + c d^2=0, a b^2 + d^3=0, a d - b c=0
- Megoldások: a=0,d=0,b c=0 -> trianguláris automorfizmus, injektív.
- Általános megoldás C felett: a^4 + b c^3=0, d= b c / a -> pl. a=1,b=1,c=-1,d=-1 -> F1=x+(x+y)^3, F2=y-(x+y)^3, automorfizmus, nincs ütközés.
- Homogén fok 3 Keller térképek osztályozása: Trace=0 és det=0 -> A3,B3 proporcionális L^3, L lineáris, c u + d v=0 -> automorfizmus.
- Lemma: Bináris homogén azonos fokú f,g-re Jac(f,g)=0 => f,g proporcionális. Bizonyítás Euler-formulával és hányados deriválttal.
- Következtetés: Homogén fok 3 Keller térképek mind automorfizmusok.

### 5. Elimináción alapuló parametrizációk és redukciók kisebb algebrai rendszerekre
Motiváció: Homogén dekompozícióval lineáris egyenletek megoldása, maradék nemlineáris top feltételek.
- J_m = (A_{m+1})_x + (B_{m+1})_y + sum_{i+j=m+2} Jac(A_i,B_j)
- J1=0 => B2 paraméterezése: b11=-2a20, b02=-a11/2, szabad: a20,a11,a02,b20 (4 dim).
- J2 lineáris A3,B3-re: b21 = -3a30 - c_x2, b12 = -a21 - c_xy/2, b03 = (-a12 - c_y2)/3, ahol c_x2,c_xy,c_y2 = Jac(A2,B2) együtthatói. Szabad: a30,a21,a12,a03,b30 (5 dim).
- J3 lineáris A4,B4-re: b31 = -4a40 - c_x3, b22 = (-3a31 - c_x2y)/2, b13 = (-2a22 - c_xy2)/3, b04 = (-a13 - c_y3)/4. Szabad: a40,a31,a22,a13,a04,b40 (6 dim).
- Top feltételek: J4: Jac(A2,B4)+Jac(A3,B3)+Jac(A4,B2)=0 (5 egyenlet), J5: Jac(A3,B4)+Jac(A4,B3)=0 (6 egyenlet), J6: Jac(A4,B4)=0 (7 egyenlet).
- J6 => B4 = k A4 vagy A4=0.
- Esetszétválasztás d=3-ra:
  * A3=0,B3≠0: Jac(A2,B3)=0 => A2=a L^2, B3=b L^3
  * B3=0,A3≠0: Jac(A3,B2)=0 => A3=a L^3, B2=b L^2
  * Mindkettő nem nulla proporcionális: A3=c L^3, B3=cp L^3, B2 = k A2 + d L^2
  * Mindkettő nulla: kvadratikus eset, nincs ütközés.
- d=3 L=x család megoldása J1,J2,ütközés alapján: a02=-1/lam^2, a11=2/lam, a20=-1, b20=-lam, c=0, cp=lam, d=0 -> A2=-(x - y/lam)^2, B2=-lam (x - y/lam)^2, B3=lam x^3, de J3 ≠0, így nem Keller.
- Általános L=x+t y család: a20=-1/(lam t+1), stb., szintén J3 ≠0.
- d=3 teljes rendszer Gröbner: 9 ismeretlen, 11 egyenlet, Groebner basis számítása folyamatban, numerikus nsolve 5 próbálkozás, maradék J6 nagy (16), nincs megoldás.
- Lemma: d=3-ra nincs megoldás, összhangban ismert tétellel, hogy harmadfokú Keller térképek 2D-ben automorfizmusok.

### 6. További családok következtetésekből
- Faktorizált ütközéses alak: A(1,0)=-1, B(1,0)=0 automatikusan teljesül, ha A(x,0)= -x^2 + x^2(x-1) v(x) + y S, B(x,0)= x^2(x-1) w(x) + y T, ahol v,w polinomok, S,T ∈ (x,y).
  * Motiváció: F1(x,0)=x(x-1) t1(x), F2(x,0)=x^2(x-1) t2(x) faktorizáció.
  * Kis fokú S,T-vel (lineáris, kvadratikus) J feltétel 7-17 egyenlet, 6-14 ismeretlen, nincs megoldás (sympy solve üres).
  * Numerikus keresés: random hill climbing 10000 iteráció, objektív 463 -> 2.9, nem nulla.
  * S,T fok 3 (9+9 coeff) + v0,v1,w0,w1 (4) =22 ismeretlen, 27 egyenlet, random search 10000 iter: objektív 3462 -> 4.8, nem nulla.
- x-fokban alacsony fokú térképek: A = p0(y)+ x p1(y)+ x^2 p2(y), B hasonló, p2(0)=-1, q2(0)=0.
  * Motiváció: x-ben lineáris/kvadratikus térképek ODE rendszerre vezetnek.
  * x-fok 1: s' (1+q) - q' s=0 => s=C(1+q), s(0)=0 => s=0, majd (1+q)(1+r')=1 => q=0, r=0 -> nincs ütközés.
  * x-fok 2, y-fok 2: 10 ismeretlen, 14 egyenlet, nincs megoldás.
  * x-fok 2, y-fok 3: 16 ismeretlen, 22 egyenlet, nsolve 10 próbálkozás, 2 közelítő megoldás de maradék -8, nem teljes.
- d=4 teljes rendszer numerikus keresés:
  * 15 ismeretlen, 20 egyenlet, nsolve 5 próbálkozás, 2 megoldás de J6 maradék 16.
  * Random hill climbing 10000 iter, objektív 6271 -> 2.99, maradék ~0.3-0.9, nem nulla.
  * A4=0 eset: 10 ismeretlen, 13 egyenlet, nsolve 10 próbálkozás, maradék 11-23, nincs megoldás.
  * Redukált A4=0 b40 eliminálva ütközéssel: 9 ismeretlen, 12 egyenlet, nsolve 20 próbálkozás, 1 megoldás de maradék 11-23.
  * Proporcionális top B4=k A4: 15 ismeretlen, 17 egyenlet, nsolve 20 próbálkozás, nincs konvergencia.
- d=5 teljes rendszer: 22 ismeretlen, 32 egyenlet, random hill climbing 5000 iter, objektív 879k -> 386, maradék több, nincs megoldás.
- d=5 ritka támogatások: nem próbálva időhiány miatt, de várhatóan nincs.

## Végrehajtott számítási módszerek
- Sympy 1.14.0: Poly, groebner, solve, nsolve, diff, Jac
- Homogén dekompozíció és lineáris elimináció kézzel levezetve
- Numerikus nsolve 5-20 random kezdőponttal, tol=1e-14, maxsteps=200-300, prec=50-80
- Random hill climbing / simulated annealing komplex együtthatókkal, 2000-10000 iteráció, lépésméret 0.5 -> 0.18 exponenciális hűtéssel
- Ritka támogatás enumeráció: kombinatorikus generálás itertools.combinations, 100-200 random minta
- Faktorizált alak paraméterezése: A=-x^2 + x^2(x-1)v(x)+ y S, B=x^2(x-1)w(x)+ y T
- x-fok szerinti osztályozás: p_i(y), q_i(y) polinomok

## Köztes lemmák (bizonyított vagy numerikusan alátámasztott)
- Lemma 1: Bináris homogén azonos fokú f,g-re Jac(f,g)=0 => f,g lineárisan összefüggő (proporcionális). Bizonyítás: Euler-formula + hányados derivált.
- Lemma 2: Homogén fok d Keller top feltétel Jac(A_d,B_d)=0 => A_d,B_d proporcionális.
- Lemma 3: Ha A_d=c L^d, B_d=cp L^d, akkor J_{2d-3}: Jac(A_{d-1},B_d)+Jac(A_d,B_{d-1})=0 => d A_{d-1} - c B_{d-1} = e L^{d-1} alakú.
- Lemma 4: d=3 esetén minden Keller térkép automorfizmus, nincs ütközés (ismert tétel, de mi is levezettük esetszétválasztással).
- Lemma 5: Faktorizált ütközéses alakban A(x,0)= -x^2 + x^2(x-1)v(x) általános alakja minden A(1,0)=-1, A∈(x,y)^2 polinomnak.
- Lemma 6: x-fokban 1 fokú Keller térképek (A=p0(y)+ x p1(y), B=q0(y)+ x q1(y)) csak triviális megoldást adnak (p1=q1=0).
- Numerikus megfigyelés: d=4,5 random keresés objektív függvénye 4000-800k-ról 3-400-ra csökken, de nem éri el a nullát, lokális minimumokba ragad, ami inkonzisztenciára utal.

## Eredmények
- d=3: teljes keresés, nincs megoldás (bizonyított esetszétválasztással és numerikusan).
- d=4: teljes (15 vars), A4=0 (10 vars), proporcionális top (15 vars), faktorizált (14-22 vars) családokban nincs megoldás 10000 iterációig.
- d=5: 22 vars, 32 egyenlet, random search 5000 iter, objektív 386, nincs megoldás.
- Ritka támogatások d=3,4: 300+ minta, 0 megoldás.
- Következtetés: d≤5-ig nincs ellenpélda a vizsgált családokban, ami összhangban van azzal, hogy a Jacobian-sejtés 2D-ben d≤4-re ismert, de nyitott magasabb fokokra.

## Nyitott irányok
- d=6,7 magasabb fokú teljes keresés 30+ változóval, 50+ egyenlettel, HPC-vel.
- Több lineáris forma L1,L2 kombinációja: A_k = sum c_{k,i} L1^{a} L2^{b}
- Gröbner-bázis számítás Singular/Magma-val nagyobb fokokra.
- Homotópia folytatás (Bertini) komplex megoldások keresésére.
- Nem-kommutatív vagy Poisson-struktúra használata.

## Új elméleti redukció: y-fok ≤2 lehetetlen (checkpoint)

A checkpoint bizonyítása szerint ha $P=p_0+p_1y+p_2y^2$, $Q=q_0+q_1y+q_2y^2$, akkor $J(P,Q)=1$ az alábbi 4 egyenletre bomlik:
- $p_0'q_1-p_1q_0'=1$
- $2p_0'q_2+p_1'q_1-p_1q_1'-2p_2q_0'=0$
- $2p_1'q_2+p_2'q_1-p_1q_2'-2p_2q_1'=0$
- $2(p_2'q_2-p_2q_2')=0$

Esetszétválasztás $p_2,q_2$ nulla/nem nulla szerint mind a 4 esetben ellentmondásra vezet a normalizációval $p_0(0)=q_0(0)=p_0(1)=q_0(1)=0$, $p_0'(0)=1$, $q_0'(0)=0$, $p_1(0)=0$, $q_1(0)=1$. Tehát minden ellenpéldában legalább az egyik komponens $y$-foka $\ge3$.

## y-fok 3 redukció

Legyen $R=Q-\lambda P=r_0+r_1y+r_2y^2$, $P=p_0+p_1y+p_2y^2+p_3y^3$, $J(P,R)=1$. Ekkor:
- $y^0: p_0'r_1-p_1r_0'=1$
- $y^1: 2p_0'r_2+p_1'r_1-p_1r_1'-2p_2r_0'=0$
- $y^2: 2p_1'r_2+p_2'r_1-p_1r_2'-2p_2r_1'-3p_3r_0'=0$
- $y^3: 2p_2'r_2+p_3'r_1-2p_2r_2'-3p_3r_1'=0$
- $y^4: 2p_3'r_2-3p_3r_2'=0\Rightarrow r_2=\alpha h^2, p_3=\beta h^3$
- Következő: $p_2=(3\beta/2\alpha)hr_1+\gamma h^2$

Maradék 3 egyenlet $p_0,p_1,r_0,r_1,h$-ra. $h$ lineáris, $r_1$ kvadratikus, $p_1$ kvadratikus, $p_0,r_0$ köbös ansatz-cal 14 ismeretlen, 19 egyenlet (2 ütközés). Numerikus keresés cmath lambdify-val 20000 iter: objektív 241k$\to$0.4, maradék 0.1-0.5, nem nulla. $h$ kvadratikus, $r_1$ köbös, $p_1$ köbös, $p_0,r_0$ kvartikus ansatz: 18 ismeretlen, 28 egyenlet, 30000 iter: 0.67 maradék.

Következtetés: $y$-fok 3 redukció kis fokú $h$-val nem ad megoldást, magasabb fokú $h$ szükséges.

## További kísérletek a checkpoint alapján

- Vonal Bezout: $p_0'q_1-p_1q_0'=1$, $\gcd(p_0',q_0')=1$, $q_0\not\equiv0$. Seed: $p_0=x-x^2$, $q_0=x^2(x-1)$, $p_1=8x$, $q_1=1+2x-12x^2$ kielégíti.
- Kiterjesztés $Q$ kvadratikus, $P$ köbös $y$-ban konstans top paraméterrel: $Q=q_0+q_1y+ay^2$, $P=p_0+p_1y+(\mu q_1+c)y^2+by^3$ $\Rightarrow$ $\mu=0,a=-4,c=-4$, majd $x^2$ obstrukció $120\neq0$.
- Fokpár szűrés: ha egy komponens foka $\le2$, automorfizmus; ha egyik fok osztja másikat, trianguláris redukció; ha $\gcd(\deg P,\deg Q)\le2$, Magnus/Applegate-Onishi. Így $\le8$ összfokig minden kizárva, első nem kizárt $(9,6)$.

## Fájlok
- search_d3.py: d=3 teljes rendszer
- test_family_d3.py: L-család
- case3_d3.py: d=3 esetszétválasztás
- druzkowski.py: Drużkowski és homogén fok 3
- search_d4_numeric.py: d=4 teljes nsolve
- search_d4_A4zero.py: d=4 A4=0
- search_d4_prop.py: B4=k A4
- factorized_ansatz.py: faktorizált kis támogatás
- sparse_enum_d3.py: ritka d=3
- x_degree_2.py: x-fok 2 család
- random_search_d4.py: random hill climbing d=4
- fast_search_d4.py: gyors d=4
- fast_search_d5.py: gyors d=5
- factorized_fast_d4.py: faktorizált fok2 S,T
- factorized_fast_d4_deg3.py: faktorizált fok3 S,T
- sparse_exact_d4.py: ritka egzakt d=4
- groebner_d4_A4zero.py: Gröbner kísérlet
- reduced_A4zero.py: redukált A4=0
- ydeg3_reduction.py: y-fok 3 redukció egzakt
- fast_ydeg3.py: y-fok 3 gyors keresés h lineáris
- fast_ydeg3_h2.py: y-fok 3 h kvadratikus
- sparse_exact_d6.py: ritka d=6
- factorized_fast_d4_deg3.py: faktorizált fok3 S,T (duplikált)
- keller_ideal Sage kód: checkpoint-ben megadott verifikátor

Minden script futtatva, eredmények logolva. Legjobb közelítő objektív $d=4$-re 2.99, $d=5$-re 386, $y$-fok 3 redukcióra 0.4-0.67, de egzakt nulla nem elérve.
