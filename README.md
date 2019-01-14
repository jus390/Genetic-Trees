# Genetic Trees
A web app that generates light adapted trees using genetic algorithms.

It optimizes trees based how many leaf poligons are directly exposed to light. The algorithm traces light rays from the light to each leaf and checks for occlusion. The final score is the number of unocluded leaves. You can also set if the leves are facing towords the light.

This repository contains two demos:
* oneTree: Scene contains one tree to optimize
* twoTrees: This scene contains two trees. It first optimizes the tree closer to the light and then the other. It shows how the other tree adats to the occlusion of the first tree.

## Running the app
To run the app simply serve each folder with your a http server.
