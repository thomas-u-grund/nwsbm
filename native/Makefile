# Native (C) Stata plugins of the nwsbm package: sbm_vem.
# Output: ../lib/plugins/<platform>/<plugin>_<platform>.plugin (the platform suffix matters for `net install`,
# which flattens all files into one folder). spi/ holds StataCorp's plugin interface (SPI 3.0, redistributed as
# published at https://www.stata.com/plugins/); it is not this package's own code.
#   make macos    universal (x86_64 + arm64) binaries, run on macOS
#   make unix     Linux (run on a Linux host or the CI runner)
# Windows (Visual Studio x64 prompt): cl /LD /O2 /DSYSTEM=STWIN32 /I spi <plugin>.c spi\stplugin.c /Fe:..\lib\plugins\windows\<plugin>_windows.plugin

CC = clang
CFLAGS = -O2 -Wall -Ispi

.PHONY: macos unix clean macos-sbm_vem unix-sbm_vem

macos: macos-sbm_vem

unix: unix-sbm_vem

macos-sbm_vem:
	mkdir -p ../lib/plugins/macos
	$(CC) -bundle -o sbm_vem.x86_64 sbm_vem.c spi/stplugin.c \
		-DSYSTEM=APPLEMAC -target x86_64-apple-macos10.11 $(CFLAGS)
	$(CC) -bundle -o sbm_vem.arm64 sbm_vem.c spi/stplugin.c \
		-DSYSTEM=APPLEMAC -target arm64-apple-macos11 $(CFLAGS)
	lipo -create -output ../lib/plugins/macos/sbm_vem_macos.plugin sbm_vem.x86_64 sbm_vem.arm64
	codesign -s - ../lib/plugins/macos/sbm_vem_macos.plugin
	rm -f sbm_vem.x86_64 sbm_vem.arm64
	@echo "Built ../lib/plugins/macos/sbm_vem_macos.plugin"

unix-sbm_vem:
	mkdir -p ../lib/plugins/unix
	gcc -shared -fPIC -DSYSTEM=OPUNIX $(CFLAGS) -pthread -o ../lib/plugins/unix/sbm_vem_unix.plugin \
		sbm_vem.c spi/stplugin.c -ldl -lpthread
	@echo "Built ../lib/plugins/unix/sbm_vem_unix.plugin"

clean:
	rm -f sbm_vem.x86_64 sbm_vem.arm64 ../lib/plugins/macos/sbm_vem_macos.plugin ../lib/plugins/unix/sbm_vem_unix.plugin
