## oppgave 3

#### a)

Netto strømmning blir større enn null ettersom vi ikke har ren diffusjon. Når potensialet er konstant er sannsynlighenene for bevegelse i positiv og negativ retning like. Dermed blir netto strømning null. Når sagtannpotensialet er aktivt derimot blir $p^+ {\not =} p^-$ for alle bindingspunkter (utenom tilfellet der $\alpha = 0.5$). Og partikklene der vil tendere mot en retning. Hvilken retning en partikkel på et gitt bindingspunkt tenderer mot vil avhenge av $\alpha$.

 Bare i tilfellet der sagtann potensialet er symmetrisk om toppunktet ($\alpha = 0.5$) vil det ikke være strømning ettersom det er symmetri i begge potensialene [skriv om]

Vi ser en større strøm i den første syklusen enn vi gjør i de senere (ca dobbelt så stor). Dette kommer av at når denne syklusen starter er alle proteinene fordelt jevn over rommet. Under starten av senere sykler derimot vil proteinene alle være bundet i en brønn i sagtann potensialet. Når de diffuserer herfra kommer færre til å forflytte seg over neste potensialtopp før sagtannpotensialet blir gjeldene.

#### b)

Strømmningen øker mot venstre med økt $T_p$ i intervallet circa $[1,500]$, etter dette minker strømningen sakte. Økningen i strømning mot venstre i det første intervallet skyldes nok at proteinene får mer tid til å diffusere, hvilket gjør at fler av dem kan stige over i neste konvekse del av potensialet. Nedgangen for $T_p > \sim 500$ skyldes nok dermed at vi når et punkt der proteinene er nært uniformt fordelt. Dermed vil påfølgende tissteg ikke bidra til å øke strømningen ved potensialbytte. Proteinene vil også bruke mer tid "fanget" i potenisalbrønnene. Dermed øker ikke netto strøming. Og ettersom $J_\text{avg} \propto \frac 1 T_p $ bidrar dette til at $J_\text{avg}$ synker.

#### c)
Den simulerte strømningen passer godt med den analytiske. Vi ser at for små $\alpha$ går $J_\text{avg}$ mot høyre. og for store $\alpha$ går strømmen mot venstre. Som nevnt i a) er strømmnigen 0 for $\alpha \approx 0.5$. Disse strømningsregningene gir mening gitt sannsynlighets funksjonene. Når $\alpha$ er liten er majoritenten av hver "sagtann" til høyre for toppen hvilket gir økt sannsynlighet for steg til høyre. Samme logikk forklarer strømning i negativ retning for $\alpha > 0.5$.


#### d)
Vi ser at for lave verdier av $\beta k $ er det dårlig sammsvar mellom mellom den simulerte og den analytiske vardien av $J_\text{avg}$. Men verdiene passer bedre jo mer $\beta k$ øker. Ettersom det i oppgave c) var tilnærmet perfekt sammsvar mellom simulert og analytisk tyder dette på at den analytiske løsningen eksisterer i grensen der den kjemiske energien i systemet er langt større enn den termiske. Altså $\beta k \gg 1$. 

#### e)
Sammenhengen mellom simulert og analytisk data plottet mot $T_p$ er relativt god. For de minste verdiene av $T_p$ ($T_p < \sim 100$)  er den simulerte strømmen mindre enn den analytiske. Utenfor dette intervallet er den simulerte strømmen større enn den analytiske. 

#### f)
Det er tydelig at den analytiske løsningen bygger på et sett med antagelser. Som nevnt i d) er det liten grad av likhet mellom modellene når den termiske energien dominerer systemet. Dermed kan en anta at den analytiske løsningen bygger på en antagelse om at systemet eksisterer i grensen $k_b T \ll K$.

Videre antar den analytiske løsningen kontinuitet, dette ser vi fra formlen hvor det integreres over kontinuerlige sannsynlighetsfordelinger. Differansen som oppstår av dette er liten og slår dessuten begge veier. Variasjonen kan minkes ved å kjøre simulasjonen med flere partikkler. 

Det er også verdt og merke seg at utrykket for $J_\text{avg}$ er forskjellig for simulert og analytisk modell. Den analytiske modellen regner et skritt til høyre som én potenisalbrønn til høyre, mens simuleringen regner hæyre skritt til å være et bindingspunkt til høyre.
Ettersom den slakere delen av potensialet er lengre enn den bratte vil simulerte proteiner ta flere steg mot bunnen av brønnen dersom de ender opp på den slake siden enn om de hadde endt opp på den bratte. Dermed blir det målte tallet for strømmen i det simulerte tilfellet større. For at disse skal bli likere kan man vekte leddene som gir positiv og negativ strømning for antall steg vi forventer at en partikkel vil ta i den retningen.