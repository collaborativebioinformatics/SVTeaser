import os
import sys
from acebinf import cmd_exe

for line in sys.stdin:
	line = line.strip()
	if os.path.exists(line + ".gz") and not "truv" in line:
		continue
	
	cmd1 = f"vcf-sort {line} | bgzip > {line}.gz && tabix {line}.gz"
	ret = cmd_exe(cmd1)
	if ret.ret_code != 0:
		print(f"error {line}")
