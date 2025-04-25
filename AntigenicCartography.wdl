task runQC {
    input{
        File titertable
        String cpu
        String mem
        String prefix
        String repsxdim
    }

    command <<<
        Rscript /usr/local/src/cartography/AC_QC.r --file ${titertable} --prefix ${prefix} --cpu ${cpu} --repsxdim ${repsxdim}
    >>>

    output{
        Array[File] qc_results = ["${prefix}_Regression.pdf","${prefix}_Dimension_Test_RMSE.csv"]
    }

    runtime {
        memory: "${mem}"
        cpu: "${cpu}"
        docker: "cameronnguyen/ac_qc:latest"
    }
}

task runAC {
    input{
        File titertable
        String cpu
        String mem
        String prefix
        String xy_lim
        String psizes
        String opacity
        String agoverprint
        String agsort
    }
    command <<<
        Rscript /usr/local/src/cartography/cartography.r \
            --input=${titertable} \
            --xy_lim=${xy_lim} \
            --prefix=${prefix} \
            --psizes=${psizes} \
            --opacity=${opacity} \
            --agoverprint=${agoverprint} \
            --agsort=${agsort} 
    >>>

    output{
        Array[File] ac_results = ["${prefix}_Cart_Distance_between_groups.csv","${prefix}_mapobj.rds","${prefix}_Cart_racmacs_coordinates.csv","${prefix}_Cart_Distance_within_ag.csv","${prefix}_Cart_Distance_within_sera.csv","${prefix}_Cartography.pdf"]
    }

    runtime {
        memory: "${mem}"
        cpu: "${cpu}"
        docker: "cameronnguyen/ac_qc:latest"
    }
}


workflow agQC {

    input {
        Array[File] titertables
        Array[String] prefixes
        String cpu
        String mem
        String repsxdim
        String xy_lim
        String psizes
        String opacity
        String agoverprint
        String agsort
    }

    scatter (index in range(length(titertables))) {
        call runAC {
            input:
                titertable = titertables[index],
                mem = mem,
                cpu = cpu,
                prefix = prefixes[index],
                xy_lim = xy_lim,
                psizes = psizes,
                opacity = opacity,
                agoverprint = agoverprint,
                agsort = agsort

        }
    }

    scatter (index in range(length(titertables))) {
        call runQC {
            input:
                titertable = titertables[index],
                mem = mem,
                cpu = cpu,
                prefix = prefixes[index],
                repsxdim= repsxdim
        }
    }

    output {
        Array[Array[File]] ac_out = runAC.ac_results
        Array[Array[File]] qc_out = runQC.qc_results
    }
}