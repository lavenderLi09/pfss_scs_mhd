# Git Workflow for `amrvac3.2_ring_average`

## Purpose

This workflow is for development in:

- `/Users/zhaoyan/Documents/codes/amrvac3.2_ring_average`

The goal is to make Ring Average development safe to iterate on, easy to roll back, and easy to compare across alternative implementations.

The core safety principle is:

- keep the original AMRVAC source tree untouched
- do all Ring Average work in the copied repository
- make small, named checkpoints that are easy to revisit

## Repository Roles

### Original Source Tree

- `/Users/zhaoyan/Documents/codes/amrvac3.2`

Rules:

- do not edit this tree for Ring Average development
- treat this as the untouched reference copy
- use it only for comparison, diffing, and fallback

### Working Repository

- `/Users/zhaoyan/Documents/codes/amrvac3.2_ring_average`

Rules:

- all Ring Average implementation work happens here
- all commits, branches, and tags for this project are made here

## Branch Model

### Stable Working Line

Primary stable branch:

- `ring-average-v1`

Meaning:

- this is the current stable working line for the project
- it should remain in a buildable, understandable state
- only merge or fast-forward work into this branch after the change is coherent enough to keep

This branch is not required to be production-ready, but it should not become a dumping ground for unrelated experiments.

### Feature / Experiment Branches

All task branches should use:

- `ra/<topic>`

Examples:

- `ra/param-plumbing`
- `ra/ring-metadata`
- `ra/cfl-only`
- `ra/cellcentered-pcm`
- `ra/cellcentered-plm`
- `ra/ppm-upgrade`
- `ra/test-off2270`
- `ra/debug-pole-instability`
- `ra/ct-exploration`

Rules:

- one branch should own one clear goal
- if a task has a different risk profile or design direction, make a new branch
- do not stack unrelated edits in the same feature branch
- if an experiment becomes messy, abandon the branch and start a fresh one from a known stable point

## Stacking Rule

We will use a simple stacked workflow.

### Default Stack Base

Unless there is a specific reason otherwise:

- start new branches from `ring-average-v1`

This keeps the branch history easy to understand.

### When Stacking on Another Feature Branch Is Allowed

You may stack a branch on top of another feature branch only when:

- the new work depends directly on unmerged work below it
- the dependency is clear and intentional
- the branch name and commit messages make the dependency obvious

Recommended stacked examples:

- `ra/cellcentered-plm` stacked on `ra/ring-metadata`
- `ra/ppm-upgrade` stacked on `ra/cellcentered-plm`

Avoid deep stacks unless they are necessary.

Preferred maximum active stack depth:

- 2 or 3 branches

If a stack becomes hard to understand, fold the lower branch into `ring-average-v1` first.

## Commit Rules

### Commit Early, But Only at Meaningful Boundaries

Create a commit when one of these is true:

- the code compiles
- a logical sub-step is finished
- a debug fix is isolated and understandable
- a test setup is added and runnable
- documentation has been updated to match the code

Do not wait until a whole large feature is complete before committing.

### Commit Message Prefixes

Use these prefixes consistently:

- `checkpoint:`
- `feat:`
- `fix:`
- `test:`
- `doc:`
- `refactor:`

Recommended meanings:

- `checkpoint:` a known-good restore point
- `feat:` new behavior or new implementation piece
- `fix:` bug fix or correction
- `test:` test case, smoke case, diagnostics, or validation tooling
- `doc:` notes, design docs, workflow docs
- `refactor:` structural cleanup without intended behavior change

Examples:

- `checkpoint: parameter plumbing compiles`
- `feat: add polar ring metadata helper`
- `feat: add cell-centered PCM ring-average pass`
- `fix: correct south-cap chunk indexing`
- `test: add off_2270 ring-average smoke case`
- `doc: add ring-average method summary`

### Checkpoint Commits

Use `checkpoint:` commits sparingly but deliberately.

A checkpoint commit should represent a state you would be happy to return to later. Typical checkpoint conditions:

- builds successfully
- or passes a smoke test
- or captures an important architectural milestone

Examples:

- `checkpoint: buildable parameter plumbing`
- `checkpoint: CFL-only prototype passes smoke`
- `checkpoint: cell-centered PCM path runs off_2270`

## Tag Rules

Tags are for high-value recovery points, not for every commit.

Use tag format:

- `ra-v1-<milestone>`

Examples:

- `ra-v1-start`
- `ra-v1-param-plumbing`
- `ra-v1-buildable`
- `ra-v1-cfl-only-smoke-pass`
- `ra-v1-cellcentered-pcm-pass`
- `ra-v1-cellcentered-plm-pass`
- `ra-v1-off2270-stable`

### When to Tag

Create a tag only if the state is worth reusing later, such as:

- the first clean baseline
- the first buildable implementation
- the first passing smoke test
- the first stable case result for a new algorithmic stage

Do not create tags for temporary debugging states.

## Push Rules

We will push regularly, but not carelessly.

### Allowed to Push

Push when one of the following is true:

- the branch has a coherent milestone
- the branch contains work worth backing up remotely
- the branch is ready for review or comparison
- the branch is part of a stack and the next branch depends on it

### Avoid Pushing

Do not push if:

- the branch is in a broken halfway state with no useful checkpoint
- the branch mixes unrelated experiments
- the commit history is too confusing to be worth sharing yet

### Push Style

Push branch-by-branch.

Recommended pattern:

1. commit a coherent milestone locally
2. verify that the branch state is understandable
3. push that branch
4. if another branch stacks on it, push the lower branch first

## Merge / Promote Rules

Promote work back into `ring-average-v1` only when:

- the branch goal is complete enough to preserve
- the code is still understandable
- there is at least basic verification for the claimed milestone

Promotion can be done by merge or by fast-forward, depending on branch shape. The important rule is not the git mechanism, but the milestone discipline.

## Recovery Rules

This workflow is designed so that rollback is always simple.

### Preferred Recovery Options

1. Switch back to a known branch
2. Check out a known checkpoint commit
3. Check out a known tag
4. Start a fresh branch from `ring-average-v1` or a stable tag

### If an Experiment Goes Bad

Do not try to rescue a deeply confused branch by piling more edits onto it.

Instead:

1. identify the last good commit or tag
2. create a new branch from there
3. re-apply only the useful parts

This is usually faster and safer than untangling a polluted branch.

## Recommended Milestone Sequence for This Project

The likely development ladder for this Ring Average work is:

1. `ra/param-plumbing`
2. `ra/ring-metadata`
3. `ra/cfl-only`
4. `ra/test-off2270`
5. `ra/cellcentered-pcm`
6. `ra/cellcentered-plm`
7. `ra/ppm-upgrade`
8. `ra/ct-exploration` only if needed later

Corresponding high-value tags are likely to be:

- `ra-v1-start`
- `ra-v1-buildable`
- `ra-v1-cfl-only-smoke-pass`
- `ra-v1-cellcentered-pcm-pass`
- `ra-v1-cellcentered-plm-pass`

## Day-to-Day Working Rules

For all later development in `amrvac3.2_ring_average`, follow these rules:

- always know which branch you are on before editing
- prefer a new `ra/<topic>` branch for risky or conceptually new work
- make small, readable commits
- create checkpoint commits at buildable or smoke-tested milestones
- tag only high-value stable states
- push coherent branch milestones, not random half-states
- do not edit the original `/Users/zhaoyan/Documents/codes/amrvac3.2`

## Practical Policy for This Project

From this point on, our working policy is:

- use this workflow for all Ring Average development in `/Users/zhaoyan/Documents/codes/amrvac3.2_ring_average`
- stack branches only when the dependency is real and explicit
- commit at meaningful checkpoints
- push branch milestones after they are understandable

This document is the default git workflow for the Ring Average experiment unless we explicitly revise it later.
