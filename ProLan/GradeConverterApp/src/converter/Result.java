package converter;

public class Result {
    private String grade;
    private String status;

    public Result(String grade, String status) {
        this.grade = grade;
        this.status = status;
    }

    public String getGrade() {
        return grade;
    }

    public String getStatus() {
        return status;
    }

    @Override
    public String toString() {
        return "Nilai Huruf: " + grade + " | Keterangan: " + status;
    }
}
