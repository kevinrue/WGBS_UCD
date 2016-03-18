# List unique folders that contain report files
find trim_galore -name '*report.txt' -exec dirname {} \; | sort | uniq
