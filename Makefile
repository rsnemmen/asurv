FC = gfortran
RM ?= rm -f
# Keep compatibility with modern gfortran's stricter argument checks.
FFLAGS ?= -std=legacy -fallow-argument-mismatch
LDFLAGS ?=
LDLIBS ?=

PROGRAM ?= asurv
SRC := asurv.f
SAMPLES_FILE := asurv.etc

.DEFAULT_GOAL := all

.PHONY: all help clean distclean test smoke smoke-gal1 smoke-gal2 smoke-gal3 \
	regression regression-gal1 regression-gal2 regression-gal3

all: $(PROGRAM)

help:
	@printf '%s\n' \
	  'Targets:' \
	  '  all        Build $(PROGRAM) (default)' \
	  '  help       Show this help text' \
	  '  clean      Remove the built executable' \
	  '  distclean  Remove the executable and common compiler scratch files' \
	  '  test       Run bundled smoke and regression tests' \
	  '  smoke      Run all bundled smoke tests' \
	  '  regression Run exact-output regression tests for bundled samples' \
	  '  smoke-gal1 Run the Kaplan-Meier sample smoke test' \
	  '  smoke-gal2 Run the two-sample sample smoke test' \
	  '  smoke-gal3 Run the bivariate sample smoke test' \
	  '  regression-gal1 Run the Kaplan-Meier exact-output regression test' \
	  '  regression-gal2 Run the two-sample exact-output regression test' \
	  '  regression-gal3 Run the bivariate exact-output regression test'

$(PROGRAM): $(SRC)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $< $(LDLIBS)

clean:
	$(RM) $(PROGRAM)

distclean: clean
	$(RM) *.o *.mod core fort.*

test: smoke regression

smoke: smoke-gal1 smoke-gal2 smoke-gal3

regression: regression-gal1 regression-gal2 regression-gal3

define run_sample_smoke
	@set -eu; \
	tmpdir=$$(mktemp -d); \
	trap 'rm -rf "$$tmpdir"' EXIT HUP INT TERM; \
	extract_section() { \
		section=$$1; \
		out=$$2; \
		awk -v section="$$section" ' \
			BEGIN { found = 0 } \
			$$0 ~ ("^\\*+[[:space:]]*" section "[[:space:]]*\\*+$$") { \
				found = 1; \
				getline; \
				while (getline > 0) { \
					if ($$0 ~ "^\\*+[[:space:]]*$$") exit; \
					print; \
				} \
			} \
			END { if (!found) exit 1 } \
		' "$(SAMPLES_FILE)" > "$$out"; \
	}; \
	extract_section "$(1).dat" "$$tmpdir/$(1).dat"; \
	extract_section "$(1).com" "$$tmpdir/$(1).com"; \
	( cd "$$tmpdir" && printf '%b' '$(2)' | "$(CURDIR)/$(PROGRAM)" >/dev/null ); \
	$(3); \
	printf '%s\n' 'smoke $(1): ok'
endef

define run_sample_regression
	@set -eu; \
	tmpdir=$$(mktemp -d); \
	trap 'rm -rf "$$tmpdir"' EXIT HUP INT TERM; \
	trim_trailing_blank_lines() { \
		infile=$$1; \
		outfile=$$2; \
		awk ' \
			{ lines[NR] = $$0 } \
			END { \
				last = NR; \
				while (last > 0 && lines[last] ~ /^[[:space:]]*$$/) last--; \
				for (i = 1; i <= last; i++) print lines[i]; \
			} \
		' "$$infile" > "$$outfile"; \
	}; \
	extract_section() { \
		section=$$1; \
		out=$$2; \
		awk -v section="$$section" ' \
			BEGIN { found = 0 } \
			$$0 ~ ("^\\*+[[:space:]]*" section "[[:space:]]*\\*+$$") { \
				found = 1; \
				getline; \
				while (getline > 0) { \
					if ($$0 ~ "^\\*+[[:space:]]*$$") exit; \
					print; \
				} \
			} \
			END { if (!found) exit 1 } \
		' "$(SAMPLES_FILE)" > "$$out"; \
	}; \
	extract_section "$(1).dat" "$$tmpdir/$(1).dat"; \
	extract_section "$(1).com" "$$tmpdir/$(1).com"; \
	extract_section "$(1).out" "$$tmpdir/$(1).expected"; \
	( cd "$$tmpdir" && printf '%b' '$(2)' | "$(CURDIR)/$(PROGRAM)" >/dev/null ); \
	trim_trailing_blank_lines "$$tmpdir/$(1).expected" "$$tmpdir/$(1).expected.trimmed"; \
	trim_trailing_blank_lines "$$tmpdir/$(1).out" "$$tmpdir/$(1).out.trimmed"; \
	diff -u "$$tmpdir/$(1).expected.trimmed" "$$tmpdir/$(1).out.trimmed"; \
	printf '%s\n' 'regression $(1): ok'
endef

smoke-gal1: $(PROGRAM)
	$(call run_sample_smoke,gal1,\n\n1\n1\ny\ngal1.com\nn\n,grep -q 'KAPLAN-MEIER ESTIMATOR' "$$tmpdir/gal1.out")

smoke-gal2: $(PROGRAM)
	$(call run_sample_smoke,gal2,\n\n1\n2\ny\ngal2.com\nn\n,grep -q 'TWO SAMPLE TEST' "$$tmpdir/gal2.out" && grep -q 'KAPLAN-MEIER ESTIMATOR' "$$tmpdir/gal2.out")

smoke-gal3: $(PROGRAM)
	$(call run_sample_smoke,gal3,\n\n2\ny\ngal3.com\nn\n,grep -q 'CORRELATION AND REGRESSION PROBLEM' "$$tmpdir/gal3.out" && grep -q 'PARAMETRIC EM ALGORITHM' "$$tmpdir/gal3.out")

regression-gal1: $(PROGRAM)
	$(call run_sample_regression,gal1,\n\n1\n1\ny\ngal1.com\nn\n)

regression-gal2: $(PROGRAM)
	$(call run_sample_regression,gal2,\n\n1\n2\ny\ngal2.com\nn\n)

regression-gal3: $(PROGRAM)
	$(call run_sample_regression,gal3,\n\n2\ny\ngal3.com\nn\n)
