Introduction
------------

The current work refers to the study of a problem wich arised during the
implementation of the IBM of Evolution and Speciation developed by
Marcus Aguiar<sup>1</sup>.

The model begins with a single species, homogeneously distributed over a
two-dimensional space, of identical individuals (genomically). This
species is composed by individuals, who find a partner reproduce
sexually, leave their offspring in the space, and die.

As the generations pass, the individuals accumulate differences, and
speciation occurs when there is no possible genetic flow between two
groups of individuals anymore.

The code is explained in detail in this Git Repository
<a href="https://github.com/IriLrnr/EvIBM" class="uri">https://github.com/IriLrnr/EvIBM</a>,
specificaly, the file DETAILED.md. In here, the results will be
analyzed.

Because the code is based on a model that has already been implemented
before, we have some results to look for:

-   After 1000 generations, the formation of 20 to 30 species.

-   The growth begins around generation 500

-   The individuals keep homogeneously distributed

In the code, to represent this population, a graph is implemented, where
the vertices corespond to individuals, and an arc exist between two dots
if the two individuals are genetically compatible, that is, if their
genome has at maximum a certain number of differences.

During the simulation, species connect and desconnect, as shown bellow
(it can be seen forward or backward)

![](../figs/report/species.png)

In the image, each set of dots of the same color compose a species. As
soon as genetic flow is stablished between a red and a yellow
individual, they become the same species.

In graph theory, a subgraph that is not connected to anyone else, is a
***maximal connected component***, as are the collection of dots of the
same color and their arcs in the image above. That is what we are going
to call a **species**. This is the genetic flow definition of
species<sup>2</sup>.

This model is spatial model, so the positioning of individuals plays a
central roll in species formation. In the current version, the species
are clustering in points in the space, and resulting in a
faster-than-expected exponential growth of the number of species - it
works as if there are physical barriers, when there is not. The result
is this behaviour of the position in the first twenty generations:

![](../gifs/clustering.gif)

While the expectations for the location is for the population to
continue homogenically distributed.

The goal of this report is to analyse graphically what is causing the
clustering, using the position gifs and the graphic of the number of
species per time. I will be adding to this report while I make tests and
find errors.

Analysis
--------

### Graphics and Gifs

The methods of the implementation of the model is described in the
[EvIBM repository](#https://github.com/IriLrnr/EvIBM). It has more than
one branch and the output is different. The branch used for this report
is the JB branch. It outputs the tables used to make the gifs and
images. Look for the git tags for the right images.

The gif for studying position is created in the `space.R` file. It
recieves the table of the form

    ##   id         x         y sp gen
    ## 1  0 63.945835 75.973488  0   0
    ## 2  1  9.348047 13.490242  0   0
    ## 3  2 52.021008  7.823214  0   0
    ## 4  3  6.990640 20.465508  0   0
    ## 5  4 46.142048 81.967728  0   0
    ## 6  5 57.331860 75.558083  0   0

where `id` is the individual, `x` and `y` are it’s coordinates, `sp` is
it’s species, and `gen` is the generation they belong, and returs a gif,
where each frame corresponds for the distribution of one generation. The
full simulation gif is

![](../gifs/complete_position_V1.gif)

It shows the movement of the individuals throughout the simulation,
indicating some problems: there is always around 1000 individuals. Why
are they clustering? Why do they only move in lines?

In the gif above, we can watch species form. How many? The R script
`species.R` recieves a table of the form

    ##   gen sp
    ## 1   0  1
    ## 2   1  1
    ## 3   2  1
    ## 4   3  1
    ## 5   4  1
    ## 6   5  1

where `gen` is the generation and `sp` is the number of species in that
generation. It outputs a simple plot

![](../figs/species/number_spp_V1.png%20=250x250)

The shape of ths graphic is right, but the curve is too steep, it is
happening before it should. The amound of species is also bigger than
the model pedicts.

Both images suggest implementation problems.

### Changes in the code

Some functions of the implementation cannot be tested outside the
simulation, and some problems, like clustering, may appear as a
“background mistake”. Those functions are the ones responsible for
reproduction, positioning, and choosing partner.

The features of the model, inside those functions, that may produce
those errors are (I think):

-   The `Reproduction` function
-   The dispersal of offspring

Or, could be the sum of various errors in the code.

So far, we have two versions of the reproduction function. They are as
follows:

    // In EvIBM, function ReproductionF(), library functions.h
    void ReproductionF (Graph G, Population progenitors, Population offspring, Parameters info)
    {   
        int focal, mate, other, i, n;

        i = 0;

        if (info->population_size < info->number_individuals) {
            for (focal = 0; focal < info->population_size; focal++) {
                if (Verify_Neighborhood (progenitors[focal]->neighborhood) < info->neighbors) {
                    mate = Choose_Mate(G, focal, progenitors, info);
                    if (mate != -1) {
                        Create_Offspring (progenitors, offspring, i, focal, mate, info);
                        i++;
                        info->population_size ++;
                    }
                }
            }
        } 

        for (focal = 0; focal < (G->U); focal++) {
            other = focal; 
            mate = -1;

            if (random_number() < 0.63 && Verify_Neighborhood (progenitors[focal]->neighborhood) > 2) {
                mate = Choose_Mate(G, focal, progenitors, info);
            }

            for (n = 0; n < 2; n++) {
                if (mate == -1) {
                    other = Choose_Mate (G, focal, progenitors, info);
                    if (other != -1)
                        mate = Choose_Mate(G, other, progenitors, info);
                }
            }

            if (mate != -1 && other != -1) {
                Create_Offspring (progenitors, offspring, i, other, mate, info);
                i++;
            }
            else {
                info->population_size --;
            }
        }
    }

This is the first one coded. It makes everyone reproduce with 67%
chance, and if not, there are two chances for a neighbor to reproduce in
place of the focal. If the population goes down, individuals with low
density can reproduce in its place.

Another alternative is

    // In EvIBM, function Reproduction(), library functions.h
    void ReproductionP (Graph G, Population progenitors, Population offspring, Parameters info)
    {   
        int focal, mate, other, i, n;
        double mu;
        unsigned int number_children; 
        int parent_population_size;

        i = 0;
        parent_population_size = info->population_size;

        mu = ((double) info->number_individuals) / (parent_population_size);

        for (focal = 0; focal < (G->U); focal++) {
            number_children = poisson(mu);
            for (n = 0; n < number_children; n++) {
                mate = Choose_Mate (G, focal, progenitors, info);
                if (mate != -1) {
                    Create_Offspring (progenitors, offspring, i, focal, mate, info);
                    i++;
                }
            }
            if (number_children > 0) 
                info->population_size += (number_children - 1);
            else {
                info->population_size --;
            }
        }
    }

This function uses the poisson distribution to sort the number of
children an individual will have, and then create them, reproducing with
different mates. It is interesting to notice that the number of 0’s are
around 36%, which matches the 67% chance of reproduction from the
`ReproductionF`. The mean of the poisson is calculated in a way that the
population always tend to `info->number_individuals`. In all the tests
bellow, the number of individuals in the simulation is 1000.

The dispersion is a characteristic of the model, each offspring has 99%
chance of being in the same spot as the focal parent, and 1% chance of
dispersal inside the radius of the parent.

    // In EvIBM, function Offspring_Position(), library functions.h
    if (random_number() <= 0.01) {
      movement_y = random_number()*info->radius;
        movement_x = random_number()*info->radius;
        if (random_number() < 0.5) {
            movement_x = -movement_x;
            movement_y = -movement_y;
        }
      //...
    }

This dispersal rate result in collapse, but increasing this rate gives
strange results. The dispersion with 1% chance of dispersal is in the
gif above. If we change it we have the distribution of the last
generation in the simulation looking like this

#### Using `ReproductionF ()`

![](../figs/position/final_distribution_V1_multi.png)

and the number of species like this.

![](../figs/species/number_spp_V1_multi.png)

The results show that something is wrong. Although patterns in ecology
are normal, the reason for this “stripes” of individuals and empty space
is at this code chunk: the dispersion should cover all the directions
around the focal. But it is clear in the code that `x` and `y` are
always both positive, or both negative, wich is just a mistake. This
will reflect on the growth of the number of species, because it changes
the number of encounters between possible partners. When dispersal is
too high, there is no formation of species, wich is predicted.

We can change the code so the dispersion could happen all around the
focal.

    // In EvIBM, function Offspring_Position(), library functions.h
    if (random_number() <= 0.01) {
      movement_y = (random_number()*2 -1)*info->radius;
      movement_x = (random_number()*2 -1)*info->radius;
      //...
    }

Changing the code, and using the “right” (1%) chance of dispersion, we
obtain the distribution below:

#### Using `ReproductionF ()`

![](../gifs/complete_position_V2.gif)

![](../figs/species/number_spp_V2.png)

The distribution is a little better, and the number of species doesn’t
change much. Changing again the disperse rate, the patters should
dissapear

![](../figs/position/final_distribution_V2_multi.png)

And they do. But it is visible that there is still some clustering.

![](../figs/species/number_spp_V2_multi.png)

Although the distribution looks better, and the number of species and
distribution are more coherent, none of the tests above gave
satisfactory enough results.

#### Using `ReproductionP ()`

I expected that changing the reproduction to poisson, fixing the mistake
in the dispersal, and running the tests, the results would be better,
but…

![](../gifs/complete_position_V3.gif)

![](../figs/species/number_spp_V3.png)

more clustering and weird behaviour of number of species. And changing
the dispersal rate does not correct the problem, because there is no
species formation:

![](../figs/position/final_distribution_V3_multi.png)

Different and more clustering patters arise. The number os species
formed, shown below, is also very different than the expected sigmoidal
distribution.

![](../figs/species/number_spp_V3_multi.png)

Looking again, carefully, into the positioning functions, I noticed a
problem.

The first function, the one that checks if an individual and a mate are
one in rage of other, in all the versions cited here, looks like this:

    // In EvIBM, function Verify_Distance(), library functions.h
    int Verify_Distance (Population progenitors, int focal, int mate, Parameters info, int increase)
        {
            int x_compatible, y_compatible, x_out_left, x_out_right, y_out_up, y_out_down;

            y_compatible = 0;
            x_compatible = 0;

            x_out_left = 0;
            x_out_right = 0;
            y_out_up = 0;
            y_out_down = 0;

            /* If an individual ratio reaches an end of the lattice, it will look on the other side, because the lattice work as a toroid */
            if (progenitors[mate]->x <= progenitors[focal]->x + info->radius + increase && progenitors[mate]->x >= progenitors[focal]->x - info->radius + increase) {
                x_compatible = 1;
            }
            if (progenitors[mate]->y <= progenitors[focal]->y + info->radius + increase && progenitors[mate]->y >= progenitors[focal]->y - info->radius + increase) {
                y_compatible = 1;
            }

            if (!x_compatible) {
                if (progenitors[focal]->x + info->radius + increase > info->lattice_width) {
                    x_out_right = progenitors[focal]->x + info->radius + increase - info->lattice_width;
                    if (progenitors[mate]->x <= x_out_right) {
                        x_compatible = 1;
                    }
                }
                else if (progenitors[focal]->x - info->radius + increase < 0) {
                    x_out_left = progenitors[focal]->x - info->radius + increase + info->lattice_width;
                    if (progenitors[mate]->x >= x_out_left) {
                        x_compatible = 1;
                    }
                }
            }

            if (!y_compatible) {
                if (progenitors[focal]->y + info->radius + increase > info->lattice_lenght) {
                    y_out_up = progenitors[focal]->y + info->radius + increase - info->lattice_lenght;
                    if (progenitors[mate]->y <= y_out_up) {
                        y_compatible = 1;
                    }
                }
                else if (progenitors[focal]->y - info->radius + increase < 0) {
                    y_out_down = progenitors[focal]->y - info->radius + increase + info->lattice_lenght;
                    if (progenitors[mate]->y >= y_out_down) {
                        y_compatible = 1;
                    }
                }
            }

            if (x_compatible && y_compatible) return 1;
            else return 0;
        }

It is looking for the partner in a square, and not in a circle, as it
should be. It has a second problem: values that should be floating
points are integers. The alternative function is the following:

    int Verify_Distance (Population progenitors, int focal, int mate, Parameters info, int increase)
        {
            double x, x0, y, y0, r;
            
            r = info->radius + increase;

            x0 = progenitors[focal]->x;
            y0 = progenitors[focal]->y;
            x = progenitors[mate]->x;
            y = progenitors[mate]->y;

            if (y0 >= info->lattice_lenght - r && y <= r)
                y = y + info->lattice_lenght;

            if (y0 <= r && y >= info->lattice_lenght - r)
                y = y - info->lattice_lenght;

            if (x0 >= info->lattice_width - r && x <= r)
                x = x + info->lattice_width;

            if (x0 <= r && x >= info->lattice_width - r)
                x = x - info->lattice_lenght;

            if ((x - x0) * (x - x0) + (y - y0) * (y - y0) <= r * r)
                return 1;
            else 
                return 0;
        }

where it is checking if the partner is inside the focal’s circle, we
obtain the following results:

#### Using `ReproductionP ()`

![](../gifs/complete_position_V3_1.gif)

![](../figs/species/number_spp_V3_1.png)

The results are the opposite of what was expected, but this is not the
only positioning function. We have the function `Offspring_Position`,
which puts the offspring in a square around the focal. It looks like
this:

    // In EvIBM, function Offspring_Position(), library functions.h
    void Offspring_Position (Population progenitors, Population offspring, int baby, int focal, Parameters info)
    {
        double movement_x, movement_y;

        movement_x = movement_y = 0;

        offspring[baby]->x = progenitors[focal]->x;
        offspring[baby]->y = progenitors[focal]->y;

        if (random_number() <= 0.01) {
            movement_y = (random_number()*2 - 1) * info->radius;
            movement_x = (random_number()*2 - 1) * info->radius;

            if (offspring[baby]->x + movement_x <= info->lattice_width && progenitors[focal]->x + movement_x >= 0)
                offspring[baby]->x += movement_x;

            else if (progenitors[focal]->x + movement_x > info->lattice_width)
                offspring[baby]->x = offspring[baby]->x + movement_x - info->lattice_width;

            else if (progenitors[focal]->x + movement_x < 0)
                offspring[baby]->x = offspring[baby]->x + movement_x + info->lattice_width;

            if (progenitors[focal]->y + movement_y <= info->lattice_lenght && progenitors[focal]->y + movement_y >= 0)
                offspring[baby]->y = offspring[baby]->y + movement_y;

            else if (progenitors[focal]->y + movement_y > info->lattice_lenght)
                offspring[baby]->y = offspring[baby]->y + movement_y - info->lattice_lenght;

            else if (progenitors[focal]->y + movement_y < 0)
                offspring[baby]->y = offspring[baby]->y + movement_y + info->lattice_lenght;
        }
    }

The way to change it is sorting the movement inside the circle, by
sorting a radius equal or less then the `info->radius`, and an angle
theta, then, extracting x and y coordinates from it.

    // In EvIBM, function Offspring_Position(), library functions.h
    void Offspring_Position (Population progenitors, Population offspring, int baby, int focal, Parameters info)
    {
        double movement_x, movement_y;
        double r, theta;

        movement_x = movement_y = 0;

        offspring[baby]->x = progenitors[focal]->x;
        offspring[baby]->y = progenitors[focal]->y;

        if (random_number() <= 0.01) {
            r = random_number() * info->radius;
            theta = rand_upto(360) + random_number();

            movement_y = sin(theta) * r;
            movement_x = cos(theta) * r;
            
      //[...]
    }

The results of this second version (V3.2, with both positioning
functions changed) looks like this:

#### Using `ReproductionP ()`

![](../gifs/complete_position_V3_2.gif)

![](../figs/species/number_spp_V3_2.png)

Clearly it is not getting better, indicating more hidden errors in the
model. The number os species should look sigmoidal, in this version it
is varying too much. Also, the concentration is worse than ever.

#### Using `ReproductionF ()`

Changing just the `Verify_Distance`, we obtain:

![](../gifs/complete_position_V2_1.gif)

![](../figs/species/number_spp_V2_1.png)

And changing also the `Offspring_Position`,

![](../gifs/complete_position_V2_2.gif)

![](../figs/species/number_spp_V2_2.png)

In this last simulation, the number of species looks sigmoidal, but it
explodes, and the positioning is still clustering, even if less, so
problem not solved.

To me, the most weird part is that the reproduction functions are very
similar in what they do, but they output very different results. The
`ReproductionP` seems like a more classy option, but tbe results are
worse in terms of curve shape and clustering, and only better in
population fluctiation and total number of species formed.

Continuing with the analysis, a problem was found in `ReproductionP`.
The population was increasing even if the baby was not actually created.
The new version looks like this:

    // In EvIBM, function ReproductionP(), library functions.h
    void ReproductionP (Graph G, Population progenitors, Population offspring, Parameters info)
        {   
            int focal, mate, other, n, baby;
            double mu;
            unsigned int number_children; 
            int parent_population_size;

            baby = 0;

            parent_population_size = info->population_size;

            mu = ((double) info->number_individuals) / (parent_population_size);

            for (focal = 0; focal < (G->U); focal++) {
                number_children = poisson (mu);
                for (n = 0; n < number_children; n++) {
                    mate = Choose_Mate (G, focal, progenitors, info);
                    if (mate != -1) {
                        Create_Offspring (progenitors, offspring, baby, focal, mate, info);
                        baby ++;
                        info->population_size ++;
                    }
                }
                info->population_size --; /* The focal dies */
            }
            printf("pop size: %d\n", info->population_size);
        }

Another problem was noticed: in the last attempt to fix the
`Offspring_Position`, the sorting of the angle was wrong, because it had
to be in radians. Changing this once more to

    // In EvIBM, function Offspring_Position(), library functions.h
    theta = random_number() * 2 * 3.14159265359;

Now, the test results:

#### Using `ReproductionP ()`

When we fix the `ReproductionP` as shown above, things do not get
better.

![](../gifs/complete_position_V3_4.gif)

![](../figs/species/number_spp_V3_4.png)

#### Using `ReproductionF ()`

Changing the just the theta angle doesn’t change much when using the
`ReproductionF`:

![](../gifs/complete_position_V2_3.gif)

![](../figs/species/number_spp_V2_3.png)

Conclusion
----------

Even if the tests did not fix my problem, they healped find problems and
fix them. It is obvious that the `ReproductionF ()` function shows
results more coherent with what was expected, even though the poisson
was a “better looking” function.

The problems found in positioning and verifying distances were solved
and they did make the program better.

### PLOT TWIST

Talking to my advisor just 2 days before the deadline for this project,
we noticed an error (in comparison to the original model). When an
individual does not reproduce (in `ReproductionF ()`), it is replaced by
a neighbor, who becomes focal. But in the original model, the position
of the offspring of this new focal is of the old focal. This preservs
the old positioning, and thoygh is a computational fix, it can be viwed
as the niche wich the offspring of the old focal would occupy.
Alternatively, because this new focal is in the range of the old focal,
it could be seen as an increase in the chance of dispersal of the
offspring.

the code fix is in `Create_Offspring ()`, it will receive both the focal
and the other (new focal).

    // In EvIBM, function Create_Offspring(), library functions.h
    void Create_Offspring (Population progenitors, Population offspring, int baby, int focal, int other, int mate, Parameters info)
    {
      int i;
        
        Offspring_Position(progenitors, offspring, baby, focal, info);

        for (i = 0; i < info->genome_size; i++) {
            if (progenitors[other]->genome[i] != progenitors[mate]->genome[i]) {
                if (rand_upto(1) == 1) {
                    offspring[baby]->genome[i] = progenitors[mate]->genome[i];
                }
                else {
                    offspring[baby]->genome[i] = progenitors[other]->genome[i];
                }
            }
            else {
                offspring[baby]->genome[i] = progenitors[mate]->genome[i];
            }
        }

        for (i = 0; i < info->genome_size; i++) {
            if (random_number() <= 0.00025) {
                mutation (offspring, baby, i);
            }
        }
    }

And the results are exactly what we expected, meaining that the rest of
the program is working as it should (there may be still minor problems).

![](../gifs/complete_position_V6.gif)

![](../figs/species/number_spp_V6.png)

Bibliography
============

1.Aguiar, M.A.M. de, Baranger, M., Baptestini, E.M., Kaufman, L. &
Bar-Yam, Y. Global patterns of speciation and diversity. *Nature*
**460**, 384–387 (2009).

2.Petit, R.J. & Excoffier, L. Gene flow and species delimitation.
*Trends in Ecology & Evolution* **24**, 386–393 (2009).
