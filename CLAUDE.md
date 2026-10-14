@AGENTS.md

# Claude Code notes
- At session start (and after any timeout): run `bash scripts/resume_status.sh` — it lists every worktree,
  the last feature-loop step committed, uncommitted work left behind, the next step, open PRs and GIAP
  Slurm jobs. Then read PROGRESS.md and continue from there.
- The main session is the lead; use the subagents in .claude/agents/ for the other roles (spec §14.1).
- You run on Stanford FarmShare inside a Slurm interactive allocation (spec §3). Source scripts/farmshare_env.sh.
- GitHub: `gh pr`, `gh issue`, `gh run` work (spec §3.10). Never push tags, never push to main.
- Commands longer than ~10 min or bigger than the allocation: sbatch from scripts/slurm/, record the job ID in PROGRESS.md, poll with squeue/sacct (spec §3.8).
- Never run work on a login node; cancel only GIAP jobs by ID; respect the R14 caps.
- Unreachable source: log it under PROGRESS.md › Blocked and continue elsewhere; never use unofficial copies (R6).
- Resumability: every feature-loop step (tests, implementation, review fixes) ends with a commit AND a
  push of its branch; subagents commit and push themselves. Parallel work lives in
  `local/worktrees/<topic>` (one branch each); in a worktree set
  `PYTHONPATH=$PWD/src:$PWD/tests/sim` and never install packages (shared venv).
- Merging: `scripts/dev/waitmerge.sh PR TITLE` squash-merges only when every CI check is present and
  passed (one snapshot); `scripts/dev/replant.sh WORKTREE FORKPOINT PR` moves a stacked branch onto main
  (retargets the PR before pushing).
- End every work unit: update PROGRESS.md, commit, push.
