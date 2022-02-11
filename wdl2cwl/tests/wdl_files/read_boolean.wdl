version 1.0
task read_boolean {
    command {
        echo true > true.txt
        echo false > false.txt
        echo True > True.txt
        echo False > False.txt
    }
    output {
        Boolean good_true = read_boolean('true.txt')
        Boolean good_false = read_boolean('false.txt')
        Boolean mixed_case_true = read_boolean('True.txt')
        Boolean mixed_case_false = read_boolean('False.txt')
    }
}

task read_bad_boolean {
    command {
        echo 1 > bad-true.txt
        echo 0 > bad-false.txt
    }
    output {
        Boolean bad_true = read_boolean('bad-true.txt')
        Boolean bad_false = read_boolean('bad-false.txt')
    }
}

task read_dynamic_boolean {
    input {
        String filename = "foobar"
    }
    command {
        echo true > ~{filename}
    }
    output {
        Boolean dynamic_true = read_boolean(filename)
    }
}
