# Mixgenes – A Genetic Algorithm for Breeding

This project contains the genetic algorithm that is used by cryptokitties to breed cats.

## Pseudo Code

```
def mixGenes(mGenes[48], sGenes[48], babyGenes[48]):
  # PARENT GENE SWAPPING
  for (i = 0; i < 12; i++):
    index = 4 * i
    for (j = 3; j > 0; j--):
      if random() < 0.25:
        swap(mGenes, index+j, index+j-1)
      if random() < 0.25:
        swap(sGenes, index+j, index+j-1)

  # BABY GENES
  for (i = 0; i < 48; i++):
    mutation = 0

    # CHECK MUTATION
    if i % 4 == 0:
      gene1 = mGene[i]
      gene2 = sGene[i]
      if gene1 > gene2:
        gene1, gene2 = gene2, gene1
      if (gene2 - gene1) == 1 and iseven(gene1):
        probability = 0.25
        if gene1 > 23:
          probability /= 2
        if random() < probability:
          mutation = (gene1 / 2) + 16

    # GIVE BABY GENES
    if mutation:
      baby[i] = mutation
    else:
      if random() < 0.5:
        babyGenes[i] = mGene[i]
      else:
        babyGenes[i] = sGene[i]
```

### Parent Gene Swapping

This is how you can get your recessive genes into the dominant slot. The best explanation is by example.

Let’s say you have a kitty with traits `(D, R1, R2, R3)`, that’s `dominant`, `recessive1`, `recessive2`, and `recessive3`. What the parent gene swapping does is with `25%` chance it will swap `R3` and `R2`, then `R2` with `R1`, then `R1` with `D`. Keep in mind that these happen one after the other, so you could end up swapping `R3` down to `R2`, then to `R1`, then to `D` all in one go. This is extremely unlikely. Like `1.5%` unlikely. Don’t expect your `R3` genes to become dominant in one step.

### Mutation

This part is probably the most difficult to follow. First, understand that this only happens after the gene swapping above.

So you may have your genes in different places by the time this happens. The first line if `i % 4 == 0` makes sure we only mutate on the dominant gene. Next we swap the genes so that the higher valued gene is `gene2`. Now we can test the precondition for a mutation if `(gene2 — gene1) == 1` and `iseven(gene1)`. Yes, that’s right, you not only need specific combinations of genes, but you also need the value to be even. Keep in mind here that `Kai` values are `+1` from their binary. All the math here should be done on the binary value.

Once the precondition passes, we do a mutation with `25%` probability if the value is less then `23`, and half of that otherwise. And the mutation is just as predictable, `mutation = (gene1 / 2) + 16`.

After that the last section is easy. Give your baby some genes! If there was a mutation, they get that. Otherwise, there’s a `50%` chance for each parent.
