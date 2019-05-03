=========================
Implementační úkol MI-MKY
=========================
Autor: Vít Souček (soucevi1@fit.cvut.cz)

Zadání
------

ECDLP v EC nad :math:`\mathbb{Z}_{p}` - Babystep-Giantstep

Mějme grupu bodů supersingulární eliptické křivky :math:`E`: :math:`y^2=x^3+1` nad tělesem
:math:`\mathbb{Z}_{p}​`, kde :math:`p=2^{42} + 1597`.
Určete diskrétní logaritmus :math:`\log_P Q`, kde :math:`P=(3,678235393584)` a :math:`Q` je dáno jako:

:math:`Q=(1528056769918,1937034018107)`

:math:`Q=(3107404331497,2558860789476)`

Implementace
------------

Úkol je implementován v jazyce Python 3.7 (pravděpodobně bude fungovat i pro Python 3.6).
Pro spuštění výpočtu stačí spustit příkaz::

   $ python3.7 ecdlp.py


Soubor ``ecdlp.py``
~~~~~~~~~~~~~~~~~~~
Soubor ``ecdlp.py`` obsahuje funkci ``main()`` a funguje jako hlavní spouštěcí modul programu.
Ve funkci ``main()`` se inicializují potřebné parametry (eliptická křivka a její body :math:`P,Q`) a následně se volá funkce ``find_logarithm`` z modulu ``bsgs``.

Pro rychlejší ladění je možné zaměnit volání funkce ``initialize()`` za funkci ``test_initialize()``,
která parametry inicializuje s menšími hodnotami a výpočet tak trvá pouze zlomek sekundy.

Soubor ``bsgs.py``
~~~~~~~~~~~~~~~~~~
Soubor ``bsgs.py`` je modul pro výpočet diskrétního logaritmu pomocí algoritmu Babystep-Giantstep.

Jeho hlavní funkci ``find_logarithm()`` je předán bod P a seznam bodů Q, jejichž logaritmus při základu P chceme najít.
Seznam bodů Q se předává najednou, aby stačilo malé kroky spočítat pouze jednou.

Funkce ``generate_baby_steps()`` vrátí seznam napočítaných malých kroků (násobků bodu P) seřazený podle souřadnic,
aby v něm šlo efektivně hledat binárním vyhledáváním. Každý z bodů :math:`X` v seznamu je vypočítán jako :math:`X=j\cdot P`
pro :math:`j \in [0, m-1]`, kde :math:`m = \sqrt{ord(P)}`.

Funkce ``giant_steps()`` počítá velké kroky a hledá, zda právě vypočítaný prvek je obsažen v seznamu malých kroků.
Velký krok je počítán jako :math:`X = Q - j \cdot m \cdot P`, kde :math:`m` je stejné jako v případě malých kroků
a :math:`j \in \mathbb N`. Po nalezení kolize platí :math:`i \cdot P = Q - j \cdot m \cdot P`
pro nějaká :math:`i,j \in \{0, ..., m-1\}`. Tedy :math:`i \cdot P +j \cdot m \cdot P = Q`,
potom :math:`(i +j \cdot m) \cdot P = Q` a :math:`\log_P Q = i+j\cdot m`.

Soubor ``elliptic_curve.py``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Tento soubor je modulem pro všechny prostředky potřebné k počítání na eliptických křivkách.

Třída ``EllipticCurve``
^^^^^^^^^^^^^^^^^^^^^^^
Třída reprezentuje jednu eliptickou křivku nad konečným tělesem :math:`F = \mathbb Z_p` ve zjednodušeném zápisu
:math:`y^2 = x^3 + ax + b`. Konstruktor bere parametry :math:`a`, :math:`b` a konečné těleso :math:`F`.

Metodou ``is_point()`` je možné ověřit, zda se bod nachází na této křivce.

Metoda ``order_approx()`` shora odhadne řád křivky podle Hasseho věty, která říká, že
:math:`p+1 -2\sqrt p \leq \#EC(GF(p)) \leq p+1 + 2\sqrt p`. Řád křivky je tedy odhadnut jako :math:`p+1+2\sqrt p`.
Se skutečným řádem křivky tento program nepočítá.

Třída ``ECPoint``
^^^^^^^^^^^^^^^^^
Třída reprezentuje jeden bod na eliptické křivce.

Konstruktor má jako parametry souřadnice :math:`x`, :math:`y` a samotnou eliptickou křivku (instanci ``EllipticCurve``).

Metoda ``order_approx()`` shora odhadne řád bodu. Řád bodu je zde shora odhadnut řádem jeho křivky (tedy též odhadem podle
Hasseho věty).

Tato třída má velké množství přetížených operátorů tak, aby šlo s bodem pohodlně počítat.
Těmito operátory jsou:

- '``+``': operace na grupě EC implementovaná podle handoutu MI-MKY
- unární '``-``': negace souřadnice :math:`y`
- binární '``-``': odčítání
- '``*``': násobení algoritmem *double-and-add*
- '``==``': rovnost párů souřadnic :math:`(x, y)` po složkách
- '``<``': nerovnost pro seřazení bodů kvůli efektivnějšímu hledání: :math:`A<B \Leftrightarrow (A.x < B.x) \lor ((A.x = B.x) \land (A.y < B.y))`

Třída ``ECPointAtInfinity``
^^^^^^^^^^^^^^^^^^^^^^^^^^^
Třída reprezentující bod v nekonečnu je podtřídou ``ECPoint``, je tedy speciálním případem bodu.
Implementována je pro pohodlnější počítání s bodem v nekonečnu.

Soubor finite_field.py
~~~~~~~~~~~~~~~~~~~~~~
Soubor je modul pokrývající potřebné operace nad konečným tělesem.

Třída ``FiniteField``
^^^^^^^^^^^^^^^^^^^^^
Třída reprezentuje jedno konečné celočíselné těleso. Nese si svůj modul a dokáže generovat svoje prvky (instance třídy
``FiniteFieldElement``).

Třída ``FiniteFieldElement``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Třída reprezentuje jeden prvek konečného tělesa. Pamatuje si svou hodnotu a číslo :math:`p` -- velikost svého konečného tělesa.

Metody této třídy jsou pouze přetížené operátory pro snadnější a přehlednější počítání.

- '``+``': Sčítání modulo :math:`p`.
- binární '``-``': Odčítání modulo :math:`p`.
- unární '``-``': Záporný prvek modulo :math:`p`.
- '``*``': Násobení modulo :math:`p`.
- '``/``': Dělení modulo :math:`p` -- Násobení inverzním prvkem vypočítaným rozšířeným Euklidovým algoritmem.
- '``**``': Mocnění modulo :math:`p`.
- '``==``': Rovnost prvků.
- '``<``': Menší než.
- '``>``': Větší než.

Soubor ``helper_tools.py``
~~~~~~~~~~~~~~~~~~~~~~~~~~
Soubor obsahující pomocné nástroje.

Třída ``BabyStepPoint`` je wrapper pro dvojici ``ECPoint`` a jeho index v seznamu malých kroků.
Index je potřeba zachovat i po seřazení, protože s na jeho základě počítá konečný výsledek logaritmu.

Funkce ``binary_search()`` binárně hledá prvek v seznamu instancí ``BabyStepPoint``.


Soubor ``vystup.txt``
~~~~~~~~~~~~~~~~~~~~~
Textový soubor obsahující výstup programu s výsledky logaritmů pro zadaný bod :math:`P` a oba body :math:`Q`.