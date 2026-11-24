# Element Climb

**A small jumping game: a lab scientist climbs a tower of lab shelves, one element per floor, and tries to get as high as possible before the rising screen catches them.**

## IN (first version)

- One HTML file, opens with a double-click, works offline.
- Keyboard only: ← → to run, Space or ↑ to jump, P to pause, R to restart.
- Our own character (a little scientist in a lab coat, drawn by the page) and our own look: lab shelves as platforms, lab-tile walls.
- Floor numbers are elements: floor 1 = H, floor 2 = He, floor 6 = C … Every 10th floor shows the full element name. After 118 it continues with numbers.
- Momentum: the longer you run, the faster you get. Faster running means higher jumps. Bouncing off a side wall keeps your speed.
- **Combo jumps:** a jump that skips 2 or more floors starts a combo. Keep landing multi-floor jumps quickly to keep it going. A short timer bar shows how much time is left. Combos give bonus points.
- **Rising pressure:** after the first few floors the screen starts scrolling up on its own and gets faster every 30 seconds. If you fall below the screen, the game is over.
- Platforms get shorter as you climb higher.
- **Score** = floors × 10 + combo bonuses. Best score is remembered in this browser.
- Start screen with the controls; game-over screen with score, floor reached, best combo and "R to play again".

## OUT (not now)

- Sound and music.
- Phone / touch controls, gamepad.
- Online high-score list.
- Characters to unlock, shop, items.

## Done when

1. You can play from the start screen to game over and restart with R.
2. Running fast gives clearly higher jumps than standing still.
3. A multi-floor jump starts a combo and the combo adds bonus points.
4. The screen scrolls by itself and gets faster over time.
5. The best score is still there after closing and reopening the page.
